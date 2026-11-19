#!/usr/bin/env python3
"""Stage, run and audit one li2025_direct CP2K/PLUMED production case.

Works for the initial 10 ps segments (bulk_r1, interface_r1), the continuation
segments (bulk_r1_seg2, interface_r1_seg2) and the fully unbiased 3 ps replicas
(event_unbiased_s1..s3 for the bulk cell, interface_event_unbiased_s1..s3 for the
15 x 15 x 30 A interface cell).

The runner is the single place that:

  * refuses an unregistered manifest or a case name that the manifest's whitelist
    entry does not list (MANIFEST_CASES below), so a typo cannot run a case the
    preparation scripts never wrote,
  * refuses a case whose declared `box_class` does not match both its
    `box_angstrom` record and the `ABC` line of the staged CP2K input
    (BOX_CLASS_WHITELIST: bulk = 14.2296 A cube, interface = 15 x 15 x 30 A),
    so an interface run can never silently inherit the bulk cell,
  * refuses a workdir that is not a direct child of --workdir-root (default the
    campaign directory), so a typo cannot write into a parent or sibling tree,
  * refuses to overwrite a workdir that already holds CP2K outputs or a
    run_record.json (stale-output gate),
  * copies the staged files (new inputs plus parent checkpoint/XYZ/PLUMED/
    HILLS/COLVAR for a continuation) under non-colliding names,
  * refuses to run when any file's SHA256 differs from the manifest, when the
    case's own .inp is not staged, or when COORD_FILE_NAME / PLUMED_INPUT_FILE
    point at a file the staging plan does not provide,
  * launches cp2k.psmp with the PLUMED shared library on LD_LIBRARY_PATH,
  * requires "PROGRAM ENDED AT" in the CP2K output,
  * for a continuation, verifies the MD step range, the trajectory range, the
    energy end step and that HILLS/COLVAR were appended rather than rewritten,
  * writes run_record.json with input hashes, step ranges, frame counts, SCF
    statistics and energies.

It does not compute a free-energy surface and must not be read as one.

Usage:
  python3 run_cp2k.py --case bulk_r1 [--ntasks 64] [--workdir DIR] [--dry-run]
  python3 run_cp2k.py --manifest seg2_manifest.json --case bulk_r1_seg2 \
      --parent-dir ../../results/li2025_direct/bulk_r1_fix1 --dry-run
  python3 run_cp2k.py --manifest interface_event_unbiased_manifest.json \
      --case interface_event_unbiased_s1 --dry-run

Exit codes: 0 success, 2 CP2K returned non-zero, 3 no PROGRAM ENDED AT,
4 staging/hash/cell/whitelist mismatch, 5 unsafe workdir path,
6 continuation or unbiased range check failed, 7 stale outputs present,
8 missing manifest or case.
"""
import argparse
import hashlib
import json
import os
import re
import shutil
import socket
import subprocess
import sys
import time
from pathlib import Path

BASE = Path(__file__).resolve().parent
CP2K = "/public/home/temp/cp2k-9.1/exe/Linux-x86-64-gfortran-mkl-plumed/cp2k.psmp"
PLUMED_LIB = "/public/home/temp/apps/plumed-2.7.3-shared/lib"
MKL_LIB = "/public/soft/intel/oneapi/mkl/latest/lib/intel64"
COMPILER_LIB = "/public/soft/intel/oneapi/compiler/latest/lib"
TBB_LIB = "/public/soft/intel/oneapi/tbb/latest/lib"
I_MPI_ROOT = "/public/soft/intel/oneapi/mpi/latest"
LEGACY_SUFFIXES = (".xyz", ".inp", ".plumed")

# Registered manifests and the only case names each may run. A manifest that is
# not listed here, or a case that its manifest does not list, is refused: the
# preparation scripts are the single source of truth for what exists.
MANIFEST_CASES = {
    "manifest.json": ("bulk_r1", "interface_r1"),
    "seg2_manifest.json": ("bulk_r1_seg2", "interface_r1_seg2"),
    "event_unbiased_manifest.json": ("event_unbiased_s1", "event_unbiased_s2",
                                     "event_unbiased_s3"),
    "interface_event_unbiased_manifest.json": ("interface_event_unbiased_s1",
                                               "interface_event_unbiased_s2",
                                               "interface_event_unbiased_s3"),
}

# Known cells, keyed by the case's `box_class`. The interface cell is NOT the
# bulk cube: a distance or minimum-image answer computed with the wrong cell is
# silently wrong, so the cell is checked against this table on every run.
BOX_CLASS_WHITELIST = {
    "bulk": (14.2296, 14.2296, 14.2296),
    "interface": (15.0, 15.0, 30.0),
}


class UnsafeWorkdirError(Exception):
    pass


class StaleOutputError(Exception):
    pass


def sha256(path):
    digest = hashlib.sha256()
    with open(str(path), "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_manifest(path=None):
    manifest_path = Path(path) if path else BASE / "manifest.json"
    if not manifest_path.is_absolute():
        manifest_path = BASE / manifest_path
    if not manifest_path.is_file():
        raise SystemExit("missing manifest: %s (run prepare_li2025_direct.py "
                         "or prepare_segment2.py)" % manifest_path)
    if manifest_path.name not in MANIFEST_CASES:
        raise SystemExit("unregistered manifest %r; register it in run_cp2k.py "
                         "MANIFEST_CASES before running it" % manifest_path.name)
    return json.loads(manifest_path.read_text(encoding="utf-8")), manifest_path.name


def get_case(manifest, name, manifest_name=None):
    allowed = MANIFEST_CASES.get(manifest_name) if manifest_name else None
    if allowed is not None and name not in allowed:
        raise SystemExit("case %r is not registered for %s (allowed: %s)"
                         % (name, manifest_name, ", ".join(allowed)))
    for record in manifest.get("cases", []):
        if record.get("name") == name:
            return record
    raise SystemExit("case %r not in manifest" % name)


def cell_lengths(text):
    """The three ABC cell lengths of a CP2K input, or None when absent."""
    match = re.search(r"^\s*ABC\s+(\S+)\s+(\S+)\s+(\S+)\s*$", text, re.MULTILINE)
    if not match:
        return None
    return [float(match.group(index)) for index in (1, 2, 3)]


def expected_box(record):
    """The whitelisted cell this case's box_class demands, or None when unset."""
    box_class = record.get("box_class")
    if box_class is None:
        return None
    if box_class not in BOX_CLASS_WHITELIST:
        raise RuntimeError("unknown box_class %r (known: %s)"
                           % (box_class, ", ".join(sorted(BOX_CLASS_WHITELIST))))
    return list(BOX_CLASS_WHITELIST[box_class])


def validate_box(record):
    """Refuse a manifest whose box_class and box_angstrom disagree."""
    expected = expected_box(record)
    if expected is None:
        return None
    declared = geometry(record)[0]
    if declared is None or len(declared) != 3:
        raise RuntimeError("case declares box_class %r but no box_angstrom"
                           % record["box_class"])
    if any(abs(float(declared[i]) - expected[i]) > 1.0e-9 for i in range(3)):
        raise RuntimeError("case declares box_class %r but box_angstrom %s != %s"
                           % (record["box_class"], list(declared), expected))
    return expected


def resolve_workdir(raw, root):
    """Require a direct child of root and reject any '..' component."""
    path = Path(raw)
    if ".." in path.parts:
        raise UnsafeWorkdirError("workdir %r contains '..'" % raw)
    workdir = path.resolve()
    root_path = Path(root).resolve()
    if workdir.parent != root_path:
        raise UnsafeWorkdirError("workdir %s is not a direct child of %s"
                                 % (workdir, root_path))
    return workdir


def output_names(case, record):
    outputs = record.get("outputs")
    project = record.get("cp2k_project")
    if outputs:
        if project and not outputs.get("trajectory", "").startswith(project + "-"):
            raise RuntimeError("case %r declares cp2k_project %r but its trajectory name is %r"
                               % (case, project, outputs.get("trajectory")))
        return outputs
    if project:
        raise RuntimeError("case %r declares cp2k_project %r but no outputs block"
                           % (case, project))
    return {
        "output": case + ".out",
        "trajectory": case + "-pos-1.xyz",
        "energy": case + "-1.ener",
        "restart": case + "-1.restart",
    }


def metad_names(case, record):
    return record.get("metad_files") or {
        "hills": "HILLS_" + case,
        "colvar": "COLVAR_" + case,
    }


def stage_entries(case, record):
    if "stage" in record:
        return record["stage"]
    return [
        {"role": "new", "source": case + suffix, "target": case + suffix,
         "sha256": record[key]}
        for suffix, key in ((".xyz", "xyz_sha256"), (".inp", "inp_sha256"),
                            (".plumed", "plumed_sha256"))
    ]


def refuse_stale(workdir, case, record):
    names = set(output_names(case, record).values())
    names.add("run_record.json")
    present = sorted(name for name in names if (workdir / name).exists())
    if present:
        raise StaleOutputError("workdir %s already contains %s; refusing to overwrite"
                               % (workdir, ", ".join(present)))


def stage_all(workdir, entries, parent_dir):
    workdir.mkdir(parents=True, exist_ok=True)
    hashes = {}
    for entry in entries:
        role = entry.get("role", "new")
        if role.startswith("parent"):
            if parent_dir is None:
                raise RuntimeError("staging %s needs --parent-dir" % entry["source"])
            source = Path(parent_dir) / entry["source"]
        elif role.startswith("new"):
            source = BASE / "inputs" / entry["source"]
        else:
            raise RuntimeError("unknown staging role %r" % role)
        if not source.is_file():
            raise RuntimeError("missing staging source: %s" % source)
        target = workdir / entry["target"]
        if not target.exists():
            shutil.copy2(str(source), str(target))
        actual = sha256(target)
        if actual != entry["sha256"]:
            raise RuntimeError("%s SHA256 mismatch: %s != %s"
                               % (entry["target"], actual, entry["sha256"]))
        hashes[entry["target"]] = actual
    return hashes


def check_case_inputs(workdir, case, record, entries):
    """The staged .inp must exist, reference only staged files, and carry the cell."""
    targets = [entry["target"] for entry in entries]
    if case + ".inp" not in targets:
        raise RuntimeError("staging plan for %r has no %s.inp" % (case, case))
    text = (workdir / (case + ".inp")).read_text(encoding="utf-8", errors="replace")
    for key in ("COORD_FILE_NAME", "PLUMED_INPUT_FILE"):
        for value in re.findall(r"^\s*%s\s+(\S+)\s*$" % key, text, re.MULTILINE):
            if Path(value).name not in targets:
                raise RuntimeError("%s references %r, which the staging plan does not provide"
                                   % (key, value))
    expected = expected_box(record)
    if expected is None:
        return None
    cell = cell_lengths(text)
    if cell is None:
        raise RuntimeError("case %r declares box_class %r but %s.inp has no ABC line"
                           % (case, record["box_class"], case))
    if any(abs(cell[i] - expected[i]) > 1.0e-9 for i in range(3)):
        raise RuntimeError("%s.inp cell %s != the %r whitelist cell %s"
                           % (case, cell, record["box_class"], expected))
    return cell


def count_frames(path):
    """Count XYZ frames written by CP2K (header line = atom count)."""
    if not path.is_file():
        return 0
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    index = 0
    frames = 0
    while index < len(lines):
        try:
            atom_count = int(lines[index].strip())
        except ValueError:
            break
        if atom_count <= 0:
            break
        index += atom_count + 2
        frames += 1
    return frames


def trajectory_steps(path):
    if not path.is_file():
        return []
    steps = []
    pattern = re.compile(r"\bi\s*=\s*(\d+)")
    with open(str(path), "r", errors="replace") as handle:
        for line in handle:
            match = pattern.search(line)
            if match:
                steps.append(int(match.group(1)))
    return steps


def energy_last_step(path):
    if not path.is_file():
        return None
    last = None
    with open(str(path), "r", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if parts and parts[0].isdigit():
                last = int(parts[0])
    return last


def table_summary(path):
    """Return (rows, last_time) for a PLUMED HILLS/COLVAR table."""
    if not path.is_file():
        return 0, None
    rows = 0
    last_time = None
    with open(str(path), "r", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            rows += 1
            try:
                last_time = float(parts[0])
            except (ValueError, IndexError):
                pass
    return rows, last_time


def parse_output(text):
    converged = [int(value) for value in
                 re.findall(r"SCF run converged in\s+(\d+)\s+steps", text)]
    steps = [int(value) for value in
             re.findall(r"MD\|\s+Step number\s+(\d+)", text)]
    energies = [float(value) for value in re.findall(
        r"ENERGY\|\s+Total FORCE_EVAL \( QS \) energy \[a\.u\.\]:?\s+(-?\d+\.\d+)", text)]
    scf = {
        "converged_runs": len(converged),
        "not_converged_runs": len(re.findall(r"SCF run NOT converged", text)),
        "iterations_min": min(converged) if converged else None,
        "iterations_max": max(converged) if converged else None,
        "iterations_mean": (sum(converged) / float(len(converged))) if converged else None,
    }
    return {
        "program_ended": "PROGRAM ENDED AT" in text,
        "abort_count": len(re.findall(r"\bABORT\b", text)),
        "md_steps_reported": len(steps),
        "md_step_first": steps[0] if steps else None,
        "md_step_last": steps[-1] if steps else None,
        "scf": scf,
        "energies_hartree": {
            "first": energies[0] if energies else None,
            "last": energies[-1] if energies else None,
            "min": min(energies) if energies else None,
            "max": max(energies) if energies else None,
        },
    }


def restart_md_start_step(path):
    """Read STEP_START_VAL from the innermost &MD block of a CP2K restart file."""
    stack = []
    values = {}
    with open(str(path), "r", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped.startswith("&"):
                parts = stripped.split()
                name = parts[0][1:].upper()
                if name == "END":
                    if stack:
                        stack.pop()
                else:
                    stack.append(name)
                continue
            if stack and stack[-1] == "MD" and stripped:
                parts = stripped.split()
                if len(parts) >= 2:
                    values.setdefault(parts[0].upper(), parts[1])
    if "STEP_START_VAL" not in values:
        raise RuntimeError("%s has no &MD STEP_START_VAL" % path)
    return int(float(values["STEP_START_VAL"]))


def geometry(record):
    geo = record.get("initial_geometry") or {}
    box = record.get("box_angstrom") or geo.get("box_angstrom")
    sc = record.get("sc_distance_angstrom") or geo.get("initial_sc_distance_angstrom")
    return box, sc


def cp2k_environment():
    env = dict(os.environ)
    parts = [PLUMED_LIB, MKL_LIB, COMPILER_LIB, TBB_LIB, env.get("LD_LIBRARY_PATH", "")]
    env["LD_LIBRARY_PATH"] = ":".join(part for part in parts if part)
    env.setdefault("I_MPI_ROOT", I_MPI_ROOT)
    env["OMP_NUM_THREADS"] = "1"
    return env


def command_for(case, ntasks, env):
    mpiexec = os.path.join(env.get("I_MPI_ROOT", I_MPI_ROOT), "bin", "mpiexec")
    return [mpiexec, "-np", str(ntasks), CP2K, "-i", case + ".inp", "-o", case + ".out"]


def audit_segment(record, parsed, frames, traj_steps, energy_step, hills, colvar,
                  parent_step):
    segment = record["segment"]
    start = segment["start_step"]
    end = segment["expected_end_step"]
    expected_first = segment["expected_first_step"]
    hills_rows, hills_last = hills
    colvar_rows, colvar_last = colvar
    baseline = record.get("metad_baseline", {})
    hills_before = baseline.get("hills_rows")
    colvar_before = baseline.get("colvar_rows")

    step_ok = (parsed["md_step_first"] == expected_first
               and parsed["md_step_last"] == end)
    traj_ok = bool(traj_steps) and traj_steps[0] > parent_step \
        and traj_steps[-1] == end
    energy_ok = energy_step == end
    overlap_ok = bool(traj_steps) and traj_steps[0] > parent_step
    hills_delta = None if hills_before is None else hills_rows - hills_before
    colvar_delta = None if colvar_before is None else colvar_rows - colvar_before
    hills_ok = hills_delta is not None and hills_delta > 0 \
        and hills_last is not None and hills_last > baseline.get("hills_last_time_fs", 0.0)
    colvar_ok = colvar_delta is not None and colvar_delta > 0 \
        and colvar_last is not None and colvar_last > baseline.get("colvar_last_time_fs", 0.0)

    return {
        "parent_case": record.get("parent_case"),
        "parent_job_id": record.get("parent_job_id"),
        "start_step": start,
        "steps": segment["steps"],
        "expected_first_step": expected_first,
        "expected_end_step": end,
        "cumulative_step_range": segment.get("cumulative_step_range"),
        "segment_time_range_fs": segment.get("segment_time_range_fs"),
        "md_step_first": parsed["md_step_first"],
        "md_step_last": parsed["md_step_last"],
        "step_check": "PASS" if step_ok else "FAIL",
        "trajectory": {
            "frames": frames,
            "first_step": traj_steps[0] if traj_steps else None,
            "last_step": traj_steps[-1] if traj_steps else None,
            "must_exceed_parent_step": parent_step,
            "expected_last_step": end,
            "check": "PASS" if traj_ok else "FAIL",
            "no_overlap_with_parent": overlap_ok,
        },
        "energy_last_step": energy_step,
        "energy_check": "PASS" if energy_ok else "FAIL",
        "hills": {
            "rows_before": hills_before,
            "rows_after": hills_rows,
            "rows_segment": hills_delta,
            "last_time_before_fs": baseline.get("hills_last_time_fs"),
            "last_time_after_fs": hills_last,
            "check": "PASS" if hills_ok else "FAIL",
        },
        "colvar": {
            "rows_before": colvar_before,
            "rows_after": colvar_rows,
            "rows_segment": colvar_delta,
            "last_time_before_fs": baseline.get("colvar_last_time_fs"),
            "last_time_after_fs": colvar_last,
            "check": "PASS" if colvar_ok else "FAIL",
        },
    }


def audit_unbiased(record, parsed, frames):
    """Step/trajectory checks for an unbiased first segment (STEPS from step 1)."""
    expected_last = int(record.get("steps", 3000))
    expected_frames = int(record.get("expected_frames", expected_last // 10 + 1))
    step_ok = parsed["md_step_first"] == 1 and parsed["md_step_last"] == expected_last
    trajectory_ok = frames == expected_frames
    return {
        "mode": "unbiased",
        "expected_first_step": 1,
        "expected_last_step": expected_last,
        "md_step_first": parsed["md_step_first"],
        "md_step_last": parsed["md_step_last"],
        "step_check": "PASS" if step_ok else "FAIL",
        "trajectory_frames": frames,
        "expected_frames": expected_frames,
        "trajectory_check": "PASS" if trajectory_ok else "FAIL",
        "velocity_policy": ("CP2K-generated MB velocities at 298.15 K from &GLOBAL SEED; "
                            "no EXT_RESTART and no inherited metaD velocities"),
        "fes": "not_applicable",
    }


def run(case, workdir, ntasks, dry_run, manifest_path, parent_dir, workdir_root):
    manifest, manifest_name = load_manifest(manifest_path)
    record = get_case(manifest, case, manifest_name)
    validate_box(record)
    workdir = resolve_workdir(workdir if workdir else str(BASE / case), workdir_root)

    if parent_dir is None and record.get("parent_workdir_name"):
        parent_dir = BASE / record["parent_workdir_name"]
    parent_dir = Path(parent_dir).resolve() if parent_dir else None
    if parent_dir is not None and parent_dir == workdir:
        raise UnsafeWorkdirError("workdir must not be the parent workdir")

    refuse_stale(workdir, case, record)
    entries = stage_entries(case, record)
    hashes = stage_all(workdir, entries, parent_dir)
    cell = check_case_inputs(workdir, case, record, entries)
    if "segment" in record:
        restart_target = next(entry["target"] for entry in entries
                              if entry.get("role") == "parent_restart")
        staged_start = restart_md_start_step(workdir / restart_target)
        if staged_start != record["segment"]["start_step"]:
            raise RuntimeError("staged restart &MD STEP_START_VAL %d != expected %d"
                               % (staged_start, record["segment"]["start_step"]))

    env = cp2k_environment()
    command = command_for(case, ntasks, env)

    if dry_run:
        payload = {
            "dry_run": True,
            "case": case,
            "manifest": manifest_name,
            "workdir": str(workdir),
            "parent_dir": str(parent_dir) if parent_dir else None,
            "command": command,
            "input_sha256": hashes,
            "staged": [entry["target"] for entry in entries],
            "box_angstrom": geometry(record)[0],
            "box_class": record.get("box_class"),
            "input_cell_angstrom": cell,
            "sc_distance_angstrom": geometry(record)[1],
        }
        if "segment" in record:
            payload["segment"] = record["segment"]
            payload["outputs"] = output_names(case, record)
            payload["metad_files"] = metad_names(case, record)
        if record.get("unbiased"):
            payload["unbiased"] = {
                "seed": record.get("seed"),
                "steps": record.get("steps"),
                "expected_frames": record.get("expected_frames"),
            }
        print(json.dumps(payload, indent=2))
        return 0

    log_path = workdir / (case + ".runner.log")
    started = time.time()
    with open(str(log_path), "w") as log:
        returncode = subprocess.call(command, cwd=str(workdir), env=env,
                                     stdout=log, stderr=subprocess.STDOUT)
    wall_seconds = time.time() - started

    outputs = output_names(case, record)
    metad = metad_names(case, record)
    output_path = workdir / outputs["output"]
    text = output_path.read_text(encoding="utf-8", errors="replace") \
        if output_path.is_file() else ""
    parsed = parse_output(text)

    trajectory = workdir / outputs["trajectory"]
    frames = count_frames(trajectory)
    traj_steps = trajectory_steps(trajectory)
    energy_step = energy_last_step(workdir / outputs["energy"])
    unbiased = bool(record.get("unbiased"))
    if unbiased:
        hills_summary = (None, None)
        colvar_summary = (None, None)
    else:
        hills_summary = table_summary(workdir / metad["hills"])
        colvar_summary = table_summary(workdir / metad["colvar"])

    result = {
        "case": case,
        "status": "COMPLETED",
        "host": socket.gethostname(),
        "job_id": os.environ.get("SLURM_JOB_ID"),
        "ntasks": ntasks,
        "cp2k": CP2K,
        "command": command,
        "workdir": str(workdir),
        "wall_seconds": round(wall_seconds, 3),
        "returncode": returncode,
        "input_sha256": hashes,
        "box_angstrom": geometry(record)[0],
        "box_class": record.get("box_class"),
        "input_cell_angstrom": cell,
        "trajectory": {
            "file": outputs["trajectory"],
            "frames": frames,
            "first_step": traj_steps[0] if traj_steps else None,
            "last_step": traj_steps[-1] if traj_steps else None,
        },
        "energy": {
            "file": outputs["energy"],
            "last_step": energy_step,
        },
        "fes": record.get("fes") or {
            "status": "not_computed",
            "note": "HILLS/COLVAR are raw metaD output; no FES is derived here.",
        },
    }
    if not unbiased:
        result["metad"] = {
            "hills_file": metad["hills"],
            "hills_rows": hills_summary[0],
            "hills_last_time_fs": hills_summary[1],
            "colvar_file": metad["colvar"],
            "colvar_rows": colvar_summary[0],
            "colvar_last_time_fs": colvar_summary[1],
        }
    if "segment" in record:
        result["segment"] = audit_segment(record, parsed, frames, traj_steps,
                                          energy_step, hills_summary, colvar_summary,
                                          record["segment"]["parent_md_step_last"])
    elif unbiased:
        result["unbiased"] = audit_unbiased(record, parsed, frames)
    else:
        result["trajectory"]["expected_frames_10ps"] = (
            record["cp2k"]["steps_initial_segment"] // record["cp2k"]["trajectory_every_fs"])
    result.update(parsed)

    failure = None
    if returncode != 0:
        failure = 2
        result["status"] = "FAILED"
    elif not parsed["program_ended"]:
        failure = 3
        result["status"] = "FAILED"
    elif "segment" in record:
        segment_audit = result["segment"]
        if segment_audit["step_check"] != "PASS" or segment_audit["energy_check"] != "PASS":
            failure = 6
            result["status"] = "FAILED_STEP_RANGE"
        elif segment_audit["trajectory"]["check"] != "PASS":
            failure = 6
            result["status"] = "FAILED_TRAJECTORY_RANGE"
        elif segment_audit["hills"]["check"] != "PASS" or segment_audit["colvar"]["check"] != "PASS":
            failure = 6
            result["status"] = "FAILED_METAD_APPEND"
    elif unbiased:
        audit = result["unbiased"]
        if audit["step_check"] != "PASS":
            failure = 6
            result["status"] = "FAILED_STEP_RANGE"
        elif audit["trajectory_check"] != "PASS":
            failure = 6
            result["status"] = "FAILED_TRAJECTORY_RANGE"

    (workdir / "run_record.json").write_text(json.dumps(result, indent=2) + "\n",
                                             encoding="utf-8")

    if failure == 2:
        print("CP2K returned %d; see %s" % (returncode, log_path))
        return 2
    if failure == 3:
        print("missing PROGRAM ENDED AT in %s" % output_path)
        return 3
    if failure == 6:
        print("continuation audit failed: %s (see %s)" % (result["status"],
                                                          workdir / "run_record.json"))
        return 6
    print("LI2025_DIRECT_OK %s frames=%d scf_not_converged=%d"
          % (case, frames, parsed["scf"]["not_converged_runs"]))
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", required=True,
                        help="case name from the manifest (bulk_r1, interface_r1, "
                             "bulk_r1_seg2, interface_r1_seg2, event_unbiased_s1..s3, "
                             "interface_event_unbiased_s1..s3); it must be registered in "
                             "MANIFEST_CASES for the chosen manifest")
    parser.add_argument("--workdir", default=None,
                        help="work directory; default <campaign>/<case>")
    parser.add_argument("--manifest", default=None,
                        help="manifest path (default manifest.json; seg2 uses "
                             "seg2_manifest.json, the unbiased replicas use "
                             "event_unbiased_manifest.json or "
                             "interface_event_unbiased_manifest.json)")
    parser.add_argument("--parent-dir", default=None,
                        help="parent workdir for continuation staging (default "
                             "<campaign>/<parent_workdir_name>)")
    parser.add_argument("--workdir-root", default=None,
                        help="the workdir must be a direct child of this directory "
                             "(default the campaign directory)")
    parser.add_argument("--ntasks", type=int,
                        default=int(os.environ.get("SLURM_NTASKS", "64")))
    parser.add_argument("--dry-run", action="store_true",
                        help="stage and verify inputs, print the command, do not run")
    args = parser.parse_args()
    root = args.workdir_root if args.workdir_root else BASE
    try:
        return run(args.case, args.workdir, args.ntasks, args.dry_run,
                   args.manifest, args.parent_dir, root)
    except UnsafeWorkdirError as error:
        print("UNSAFE WORKDIR: %s" % error)
        return 5
    except StaleOutputError as error:
        print("STALE OUTPUT: %s" % error)
        return 7
    except RuntimeError as error:
        print("STAGING/HASH ERROR: %s" % error)
        return 4


if __name__ == "__main__":
    sys.exit(main())

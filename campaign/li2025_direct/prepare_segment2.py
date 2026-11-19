#!/usr/bin/env python3
"""Build the second 10 ps segment (seg2) of the Li 2025 direct-addition replicas.

The two first segments finished normally and were pulled back:

  bulk_r1      job 8562691, 10000 steps -> results/li2025_direct/bulk_r1_fix1
  interface_r1 job 8562695, 10000 steps -> results/li2025_direct/interface_r1_fix1

This script does NOT recompute or re-run anything. It validates the terminal
parent state, then emits, inside campaign/li2025_direct/:

  inputs/<case>_seg2.inp      parent CP2K input + one inserted &EXT_RESTART block
  inputs/<case>_seg2.plumed   parent PLUMED input + one inserted RESTART action
  seg2_manifest.json          parent/new hashes, counters, staging plan, scope

Continuation semantics (all hard-gated, nothing guessed):

  * CP2K 9.1: &EXT_RESTART is a top-level section placed immediately after
    &END GLOBAL, with RESTART_DEFAULT .FALSE. and only POS/VEL/COUNTERS on.
  * STEPS stays 10000. CP2K's STEPS is the per-segment budget, so a segment
    starting at STEP_START_VAL=10000 runs steps 10001..20000. This matches the
    project's own prior art (steered_ext_final/STEERED_EXTENSION.json: parent
    start 3000 + child STEPS 1000 -> last step 4000) and metad_resume/seg1.
  * PLUMED: RESTART is inserted as the first action after UNITS (PLUMED requires
    UNITS to be its first action). HILLS/COLVAR keep their parent file names and
    are appended, never restarted from zero.
  * The parent checkpoint has TWO STEP_START_VAL values: &MD 10000 and
    &METADYN 0. Only the &MD value is authoritative; the parser is scoped to
    &MD (the project already hit this trap at restart-chain boundaries).

All parent files are staged into the new workdir under new names (the restart
becomes <case>_parent.restart), so nothing in the parent directory is modified
and no output name collides with the copied checkpoint.

Run: python3 prepare_segment2.py [--parent-root DIR] [--only CASE]
"""
import argparse
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"
DEFAULT_PARENT_ROOT = ROOT / "results" / "li2025_direct"
REMOTE_ROOT = "/public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/li2025_direct"

SEGMENT_STEPS = 10000
EXPECTED_START_STEP = 10000
EXPECTED_END_STEP = EXPECTED_START_STEP + SEGMENT_STEPS
ATOM_COUNT = 296

CASES = [
    {
        "name": "bulk_r1_seg2",
        "parent_case": "bulk_r1",
        "parent_workdir": "bulk_r1_fix1",
        "parent_job": "8562691",
        "project": "bulk_r1",
        "hills": "HILLS_bulk_r1",
        "colvar": "COLVAR_bulk_r1",
    },
    {
        "name": "interface_r1_seg2",
        "parent_case": "interface_r1",
        "parent_workdir": "interface_r1_fix1",
        "parent_job": "8562695",
        "project": "interface_r1",
        "hills": "HILLS_interface_r1",
        "colvar": "COLVAR_interface_r1",
    },
]

EXT_RESTART_TEMPLATE = (
    "&EXT_RESTART\n"
    "  RESTART_FILE_NAME {restart}\n"
    "  RESTART_DEFAULT .FALSE.\n"
    "  RESTART_POS .TRUE.\n"
    "  RESTART_VEL .TRUE.\n"
    "  RESTART_COUNTERS .TRUE.\n"
    "&END EXT_RESTART\n"
    "\n"
)


def sha256(path):
    digest = hashlib.sha256()
    with open(str(path), "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_text(path):
    return path.read_text(encoding="utf-8")


def md_block_counters(text):
    """Return the key/value pairs of the innermost &MD section.

    The restart file also contains a &METADYN block with STEP_START_VAL 0;
    scoping to the innermost section keeps that decoy out of the gate.
    """
    stack = []
    values = {}
    for line in text.splitlines():
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
        if stack and stack[-1] == "MD" and stripped and not stripped.startswith("#"):
            parts = stripped.split()
            if len(parts) >= 2:
                values.setdefault(parts[0].upper(), parts[1])
    return values


def all_values(text, key):
    pattern = re.compile(r"^\s*" + re.escape(key) + r"\s+(\S+)", re.MULTILINE)
    return pattern.findall(text)


def section_atom_lines(text, section):
    """Count atom records in a restart &COORD (element + xyz) or &VELOCITY (xyz) block."""
    stack = []
    count = 0
    section = section.upper()
    for line in text.splitlines():
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
        if not stack or stack[-1] != section:
            continue
        parts = stripped.split()
        if section == "COORD":
            if len(parts) == 4 and parts[0].isalpha():
                count += 1
        elif len(parts) == 3:
            try:
                [float(value) for value in parts]
                count += 1
            except ValueError:
                pass
    return count


def table_summary(path):
    rows = 0
    last_time = None
    for line in read_text(path).splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if not parts:
            continue
        rows += 1
        try:
            last_time = float(parts[0])
        except ValueError:
            pass
    return rows, last_time


def energy_last_step(path):
    last = None
    for line in read_text(path).splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if parts and parts[0].isdigit():
            last = int(parts[0])
    return last


def trajectory_steps(path):
    steps = []
    pattern = re.compile(r"\bi\s*=\s*(\d+)")
    with open(str(path), "r") as handle:
        for line in handle:
            match = pattern.search(line)
            if match:
                steps.append(int(match.group(1)))
    return steps


def steps_line_value(inp_text):
    matches = re.findall(r"^\s*STEPS\s+(\d+)\s*$", inp_text, re.MULTILINE)
    if len(matches) != 1:
        raise SystemExit("parent input must contain exactly one 'STEPS <n>' line, found %d"
                         % len(matches))
    return int(matches[0])


def insert_ext_restart(inp_text, restart_name):
    """Insert the &EXT_RESTART block immediately after &END GLOBAL."""
    if "/" in restart_name or "\\" in restart_name or ".." in restart_name:
        raise ValueError("restart file name must be a bare basename: %r" % restart_name)
    if "&EXT_RESTART" in inp_text.upper():
        raise ValueError("input already contains &EXT_RESTART; refusing to insert a second one")
    lines = inp_text.splitlines(keepends=True)
    ends = [index for index, line in enumerate(lines)
            if line.strip().upper() == "&END GLOBAL"]
    if len(ends) != 1:
        raise ValueError("input must contain exactly one top-level '&END GLOBAL', found %d"
                         % len(ends))
    block = EXT_RESTART_TEMPLATE.format(restart=restart_name)
    lines.insert(ends[0] + 1, block)
    text = "".join(lines)
    if text.replace(block, "") != inp_text:
        raise ValueError("insertion changed the parent input beyond the &EXT_RESTART block")
    if text.count("&EXT_RESTART") != 1 or text.count("&END EXT_RESTART") != 1:
        raise ValueError("inserted &EXT_RESTART block is not unique")
    return text


def insert_plumed_restart(plumed_text):
    """Insert RESTART as the first action after UNITS (or as line 1 without UNITS)."""
    lines = plumed_text.splitlines(keepends=True)
    if any(line.strip().upper() == "RESTART" for line in lines):
        return plumed_text, False
    units = [index for index, line in enumerate(lines)
             if line.strip().upper().startswith("UNITS ")]
    position = units[0] + 1 if units else 0
    lines.insert(position, "RESTART\n")
    return "".join(lines), True


def resolve_parent_root(argument):
    return Path(argument).resolve() if argument else DEFAULT_PARENT_ROOT


def initial_geometry(parent_case):
    """Copy the frozen box / initial S-C distance from the r1 manifest, if present."""
    path = HERE / "manifest.json"
    if not path.is_file():
        return {}
    manifest = json.loads(read_text(path))
    for record in manifest.get("cases", []):
        if record.get("name") == parent_case:
            return {
                "box_angstrom": record.get("box_angstrom"),
                "initial_sc_distance_angstrom": record.get("sc_distance_angstrom"),
            }
    return {}


def _relative_or_absolute(path):
    try:
        return str(Path(path).resolve().relative_to(ROOT))
    except ValueError:
        return str(Path(path).resolve())


def inspect_parent(case, parent_root):
    parent_dir = parent_root / case["parent_workdir"]
    if not parent_dir.is_dir():
        raise SystemExit("missing parent workdir: %s" % parent_dir)

    paths = {
        "run_record": parent_dir / "run_record.json",
        "restart": parent_dir / ("%s-1.restart" % case["parent_case"]),
        "inp": parent_dir / ("%s.inp" % case["parent_case"]),
        "plumed": parent_dir / ("%s.plumed" % case["parent_case"]),
        "xyz": parent_dir / ("%s.xyz" % case["parent_case"]),
        "hills": parent_dir / case["hills"],
        "colvar": parent_dir / case["colvar"],
        "energy": parent_dir / ("%s-1.ener" % case["parent_case"]),
        "trajectory": parent_dir / ("%s-pos-1.xyz" % case["parent_case"]),
    }
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise SystemExit("missing parent files: %s" % ", ".join(missing))

    record = json.loads(read_text(paths["run_record"]))
    if record.get("status") != "COMPLETED":
        raise SystemExit("%s parent run_record status is %r, expected COMPLETED"
                         % (case["name"], record.get("status")))
    if record.get("program_ended") is not True:
        raise SystemExit("%s parent run_record program_ended is not true" % case["name"])
    if record.get("returncode") != 0:
        raise SystemExit("%s parent returncode is %r" % (case["name"], record.get("returncode")))
    if str(record.get("job_id")) != case["parent_job"]:
        raise SystemExit("%s parent job id %r != expected %s"
                         % (case["name"], record.get("job_id"), case["parent_job"]))
    if record.get("md_step_last") != EXPECTED_START_STEP:
        raise SystemExit("%s parent md_step_last %r != %d"
                         % (case["name"], record.get("md_step_last"), EXPECTED_START_STEP))

    for name, digest in record.get("input_sha256", {}).items():
        actual = sha256(parent_dir / name)
        if actual != digest:
            raise SystemExit("%s parent input %s sha mismatch: %s != %s"
                             % (case["name"], name, actual, digest))

    inp_text = read_text(paths["inp"])
    plumed_text = read_text(paths["plumed"])
    restart_text = read_text(paths["restart"])

    if not re.search(r"^\s*PROJECT\s+%s\s*$" % re.escape(case["project"]),
                     inp_text, re.MULTILINE):
        raise SystemExit("%s parent input PROJECT does not name %s"
                         % (case["name"], case["project"]))
    if "PLUMED_INPUT_FILE ./%s.plumed" % case["parent_case"] not in inp_text:
        raise SystemExit("%s parent input does not point at ./%s.plumed"
                         % (case["name"], case["parent_case"]))
    if steps_line_value(inp_text) != SEGMENT_STEPS:
        raise SystemExit("%s parent input STEPS is not %d" % (case["name"], SEGMENT_STEPS))
    counters = md_block_counters(restart_text)
    if "STEP_START_VAL" not in counters:
        raise SystemExit("%s restart has no &MD STEP_START_VAL" % case["name"])
    start_step = int(float(counters["STEP_START_VAL"]))
    if start_step != EXPECTED_START_STEP:
        raise SystemExit("%s restart &MD STEP_START_VAL %d != %d"
                         % (case["name"], start_step, EXPECTED_START_STEP))
    if "TIME_START_VAL" not in counters:
        raise SystemExit("%s restart has no &MD TIME_START_VAL" % case["name"])
    time_start = float(counters["TIME_START_VAL"])
    if abs(time_start - EXPECTED_START_STEP) > 1.0e-3:
        raise SystemExit("%s restart &MD TIME_START_VAL %.6f is not %d fs"
                         % (case["name"], time_start, EXPECTED_START_STEP))
    restart_steps = int(float(counters.get("STEPS", "0")))
    if restart_steps != SEGMENT_STEPS:
        raise SystemExit("%s restart &MD STEPS %d != %d"
                         % (case["name"], restart_steps, SEGMENT_STEPS))

    coord_atoms = section_atom_lines(restart_text, "COORD")
    velocity_atoms = section_atom_lines(restart_text, "VELOCITY")
    if coord_atoms != ATOM_COUNT:
        raise SystemExit("%s restart &COORD has %d atom lines, expected %d"
                         % (case["name"], coord_atoms, ATOM_COUNT))
    if velocity_atoms != ATOM_COUNT:
        raise SystemExit("%s restart &VELOCITY has %d atom lines, expected %d"
                         % (case["name"], velocity_atoms, ATOM_COUNT))

    all_step_values = all_values(restart_text, "STEP_START_VAL")
    hills_rows, hills_last_time = table_summary(paths["hills"])
    colvar_rows, colvar_last_time = table_summary(paths["colvar"])
    energy_step = energy_last_step(paths["energy"])
    traj_steps = trajectory_steps(paths["trajectory"])
    if hills_rows <= 0 or hills_last_time is None:
        raise SystemExit("%s parent HILLS is empty" % case["name"])
    if colvar_rows <= 0 or colvar_last_time is None:
        raise SystemExit("%s parent COLVAR is empty" % case["name"])
    if energy_step != EXPECTED_START_STEP:
        raise SystemExit("%s parent energy ends at step %r, expected %d"
                         % (case["name"], energy_step, EXPECTED_START_STEP))
    if not traj_steps or traj_steps[-1] != EXPECTED_START_STEP:
        raise SystemExit("%s parent trajectory ends at step %r, expected %d"
                         % (case["name"], traj_steps[-1] if traj_steps else None,
                            EXPECTED_START_STEP))
    if hills_last_time >= EXPECTED_START_STEP:
        raise SystemExit("%s parent last hill time %.3f is not before %d"
                         % (case["name"], hills_last_time, EXPECTED_START_STEP))

    return {
        "parent_dir": parent_dir,
        "paths": paths,
        "record": record,
        "inp": inp_text,
        "plumed": plumed_text,
        "restart": restart_text,
        "counters": counters,
        "start_step": start_step,
        "time_start": time_start,
        "restart_steps": restart_steps,
        "coord_atoms": coord_atoms,
        "velocity_atoms": velocity_atoms,
        "all_step_values": all_step_values,
        "hills_rows": hills_rows,
        "hills_last_time": hills_last_time,
        "colvar_rows": colvar_rows,
        "colvar_last_time": colvar_last_time,
        "energy_step": energy_step,
        "traj_frames": len(traj_steps),
        "traj_last_step": traj_steps[-1],
    }


def build_case(case, parent, parent_root):
    parent_case = case["parent_case"]
    restart_target = "%s_parent.restart" % case["name"]

    new_inp = insert_ext_restart(parent["inp"], restart_target)
    new_plumed, restart_added = insert_plumed_restart(parent["plumed"])
    if steps_line_value(new_inp) != SEGMENT_STEPS:
        raise SystemExit("%s: generated input changed STEPS" % case["name"])

    inp_path = INPUTS / ("%s.inp" % case["name"])
    plumed_path = INPUTS / ("%s.plumed" % case["name"])
    inp_path.write_text(new_inp, encoding="utf-8")
    plumed_path.write_text(new_plumed, encoding="utf-8")

    stage = [
        {"role": "new_input", "source": inp_path.name, "target": inp_path.name},
        {"role": "new_plumed", "source": plumed_path.name, "target": "%s.plumed" % parent_case},
        {"role": "parent_restart", "source": parent["paths"]["restart"].name,
         "target": restart_target},
        {"role": "parent_xyz", "source": parent["paths"]["xyz"].name,
         "target": parent["paths"]["xyz"].name},
        {"role": "parent_input", "source": parent["paths"]["inp"].name,
         "target": "parent_%s.inp" % parent_case},
        {"role": "parent_plumed", "source": parent["paths"]["plumed"].name,
         "target": "parent_%s.plumed" % parent_case},
        {"role": "parent_run_record", "source": parent["paths"]["run_record"].name,
         "target": "parent_run_record.json"},
        {"role": "parent_hills", "source": parent["paths"]["hills"].name,
         "target": parent["paths"]["hills"].name},
        {"role": "parent_colvar", "source": parent["paths"]["colvar"].name,
         "target": parent["paths"]["colvar"].name},
    ]
    for entry in stage:
        if entry["role"].startswith("parent"):
            source_path = parent["parent_dir"] / entry["source"]
        else:
            source_path = INPUTS / entry["source"]
        entry["sha256"] = sha256(source_path)
        entry["bytes"] = source_path.stat().st_size

    parent_evidence = {
        "run_record": {
            "path": str(parent["paths"]["run_record"]),
            "sha256": sha256(parent["paths"]["run_record"]),
            "status": parent["record"]["status"],
            "program_ended": parent["record"]["program_ended"],
            "job_id": parent["record"]["job_id"],
            "md_step_last": parent["record"]["md_step_last"],
            "trajectory_frames": parent["record"]["trajectory"]["frames"],
            "scf_not_converged_runs": parent["record"]["scf"]["not_converged_runs"],
            "energy_last_hartree": parent["record"]["energies_hartree"]["last"],
        },
        "restart": {
            "path": str(parent["paths"]["restart"]),
            "sha256": sha256(parent["paths"]["restart"]),
            "bytes": parent["paths"]["restart"].stat().st_size,
            "md_step_start_val": parent["start_step"],
            "md_time_start_val_fs": parent["time_start"],
            "md_steps_field": parent["restart_steps"],
            "coord_atom_lines": parent["coord_atoms"],
            "velocity_atom_lines": parent["velocity_atoms"],
            "all_step_start_val_values": parent["all_step_values"],
            "note": ("the &METADYN copy STEP_START_VAL 0 is not authoritative; "
                     "only the &MD block counter drives the segment"),
        },
        "inputs": {
            "inp_sha256": sha256(parent["paths"]["inp"]),
            "plumed_sha256": sha256(parent["paths"]["plumed"]),
            "xyz_sha256": sha256(parent["paths"]["xyz"]),
        },
        "trajectory": {
            "path": str(parent["paths"]["trajectory"]),
            "sha256": sha256(parent["paths"]["trajectory"]),
            "frames": parent["traj_frames"],
            "first_step": 0,
            "last_step": parent["traj_last_step"],
        },
        "energy": {
            "path": str(parent["paths"]["energy"]),
            "sha256": sha256(parent["paths"]["energy"]),
            "last_step": parent["energy_step"],
        },
        "metad": {
            "hills_file": parent["paths"]["hills"].name,
            "hills_rows": parent["hills_rows"],
            "hills_last_time_fs": parent["hills_last_time"],
            "colvar_file": parent["paths"]["colvar"].name,
            "colvar_rows": parent["colvar_rows"],
            "colvar_last_time_fs": parent["colvar_last_time"],
        },
    }

    return {
        "name": case["name"],
        "parent_case": parent_case,
        "parent_job_id": case["parent_job"],
        "parent_workdir_name": case["parent_workdir"],
        "parent_workdir_remote": "%s/%s" % (REMOTE_ROOT, case["parent_workdir"]),
        "parent_workdir_local": _relative_or_absolute(parent["parent_dir"]),
        "cp2k_project": case["project"],
        "initial_geometry": initial_geometry(parent_case),
        "segment": {
            "start_step": parent["start_step"],
            "steps": SEGMENT_STEPS,
            "expected_first_step": parent["start_step"] + 1,
            "expected_end_step": EXPECTED_END_STEP,
            "segment_step_range": [parent["start_step"] + 1, EXPECTED_END_STEP],
            "cumulative_step_range": [1, EXPECTED_END_STEP],
            "segment_time_range_fs": [parent["time_start"], parent["time_start"] + SEGMENT_STEPS],
            "parent_md_step_last": parent["record"]["md_step_last"],
        },
        "stage": stage,
        "outputs": {
            "output": "%s.out" % case["name"],
            "trajectory": "%s-pos-1.xyz" % case["project"],
            "energy": "%s-1.ener" % case["project"],
            "restart": "%s-1.restart" % case["project"],
        },
        "metad_files": {
            "hills": parent["paths"]["hills"].name,
            "colvar": parent["paths"]["colvar"].name,
        },
        "plumed_restart_added": restart_added,
        "metad_baseline": {
            "hills_rows": parent["hills_rows"],
            "hills_last_time_fs": parent["hills_last_time"],
            "colvar_rows": parent["colvar_rows"],
            "colvar_last_time_fs": parent["colvar_last_time"],
        },
        "parent_evidence": parent_evidence,
        "fes": {
            "status": "not_computed",
            "note": ("HILLS/COLVAR are appended raw metaD output; 20 ps over two "
                     "single replicas is still not a converged FES."),
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent-root", default=None,
                        help="local root holding the pulled parent workdirs "
                             "(default results/li2025_direct)")
    parser.add_argument("--only", default=None,
                        help="build only this case name (bulk_r1_seg2 or interface_r1_seg2)")
    args = parser.parse_args()

    parent_root = resolve_parent_root(args.parent_root)
    selected = [case for case in CASES if args.only in (None, case["name"])]
    if not selected:
        raise SystemExit("unknown --only %r" % args.only)

    INPUTS.mkdir(parents=True, exist_ok=True)
    built = []
    for case in selected:
        parent = inspect_parent(case, parent_root)
        built.append(build_case(case, parent, parent_root))

    manifest = {
        "schema": "li2025_direct_seg2/v1",
        "session": "li2025-direct-pathway3-seg2",
        "purpose": ("Second 10 ps segment of the Li 2025 pathway-3 direct C-S addition, "
                    "continued from the completed first-segment checkpoints."),
        "scope": {
            "what_this_is": [
                "Continuation segments 10001-20000 (10 ps) for bulk_r1 and interface_r1.",
                "CP2K input is the parent input plus one &EXT_RESTART block; every "
                "scientific setting (functional, box, charge, thermostat, timestep, "
                "output cadence, STEPS=10000) is byte-identical.",
                "PLUMED input is the parent input plus one RESTART action; same CV, "
                "SIGMA, PACE, HEIGHT, wall and HILLS/COLVAR file names.",
            ],
            "what_this_is_not": [
                "NOT a rerun of the first segment: parent files are copied, never modified.",
                "NOT Li's three-replica statistics; still one bulk and one interface replica.",
                "NOT a converged FES: 20 ps of single-replica metaD is not a free-energy surface.",
                "NOT the +2 K+ neutral model: the q = -2, counterion-free periodic model "
                "is unchanged from the first segment.",
            ],
        },
        "continuation_policy": {
            "cp2k": ("&EXT_RESTART immediately after &END GLOBAL; RESTART_DEFAULT .FALSE.; "
                     "RESTART_POS/VEL/COUNTERS .TRUE.; RESTART_FILE_NAME points at the "
                     "copied checkpoint basename inside the segment workdir"),
            "steps_semantics": ("STEPS=10000 is the per-segment budget. A segment starting at "
                                "STEP_START_VAL=10000 runs steps 10001-20000. Verified against "
                                "steered_ext_final/STEERED_EXTENSION.json (start 3000 + STEPS "
                                "1000 -> last step 4000) and metad_resume/seg1."),
            "plumed": ("RESTART is the first action after UNITS; PLUMED appends to the staged "
                       "HILLS/COLVAR instead of restarting hill deposition from zero"),
            "staging": ("the runner copies parent artifacts into the segment workdir under "
                        "non-colliding names and verifies every SHA256 before launching CP2K"),
            "stale_outputs": ("a workdir containing a CP2K output, trajectory, energy, restart "
                              "or run_record.json is refused instead of overwritten"),
        },
        "remote_root": REMOTE_ROOT,
        "parent_root_local": str(parent_root.relative_to(ROOT)),
        "cases": built,
    }
    (HERE / "seg2_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n",
                                             encoding="utf-8")

    summary = {
        "cases": [
            {
                "name": entry["name"],
                "parent_job": entry["parent_job_id"],
                "start_step": entry["segment"]["start_step"],
                "segment_steps": entry["segment"]["steps"],
                "expected_end_step": entry["segment"]["expected_end_step"],
                "restart": entry["parent_evidence"]["restart"]["path"],
            }
            for entry in built
        ]
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Self-test for the unbiased Li 2025 INTERFACE event replicas.

Uses the real pulled interface trajectories (both 10 ps segments) for
source/integration checks and synthetic files for the hard gates. No cluster,
no CP2K, no HPC submission.

Covered:
  * frame identification by MD step + numeric time inside 0.01 fs: the XYZ comment
    prints 10410.000 while the COLVAR/report value is 10409.999160, so a
    string-precision assumption would fail; zero or several matches must stop,
  * verbatim atom lines, composition, atom order, intact waters/HCHO/SO3, contact
    gate and S-C/C=O/C-H/S-O distances computed in the 15 x 15 x 30 A cell,
  * an independent re-scan of all 2001 frames: exactly six frames below 2.0 A,
    one contiguous window (steps 10370-10420), shortest at step 10410,
  * COLVAR cross-check at that step,
  * unbiased input: no PLUMED/METADYN/FREE_ENERGY/EXT_RESTART/VELOCITY, STEPS 3000,
    ABC 15.0 15.0 30.0 and no 14.2296 anywhere, three replicas differing only in
    &GLOBAL SEED,
  * six wrappers (compute + gpu partition): 64 ranks, 32G, 48 h, mail END,FAIL,
    no --gres, distinct workdirs/job names/logs, analyzer run with --box 15,15,30,
  * runner: registered manifest/case whitelist, box_class vs box_angstrom and vs
    the staged ABC line, staged COORD/PLUMED references, dry-run staging, hash
    mismatch (exit 4), stale output (exit 7), unsafe workdir (exit 5), complete
    and truncated fake runs (exit 0 / 6),
  * analyzer: --box is mandatory, the interface cell reproduces S-C 1.808134 A,
    no verdict, TSV output,
  * Python 3.6 syntax/API and `bash -n` for every script in the directory.

Run: python3 selftest_interface_event_unbiased.py
Exit 0 and a final SELFTEST_INTERFACE_EVENT_UNBIASED_OK line on success.
"""
import sys

sys.dont_write_bytecode = True

import ast
import difflib
import hashlib
import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))

import prepare_event_unbiased as bulk
import prepare_interface_event_unbiased as prep
import run_cp2k as runner

MANIFEST = HERE / "interface_event_unbiased_manifest.json"
MANIFEST_NAME = "interface_event_unbiased_manifest.json"
INTERFACE_BOX = "15.0,15.0,30.0"
CHECKS = {"count": 0}


def ok(label, condition, detail=""):
    CHECKS["count"] += 1
    if not condition:
        raise AssertionError("%s%s" % (label, (": " + detail) if detail else ""))
    print("ok %2d - %s" % (CHECKS["count"], label))


def read_text(path):
    return Path(path).read_text(encoding="utf-8")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run_process(args):
    process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = process.communicate()[0].decode("utf-8", "replace")
    return process.returncode, output


def run_runner(args):
    return run_process([sys.executable, str(HERE / "run_cp2k.py")] + args)


def independent_sc_scan(paths, box):
    """Independent minimum-image d(S-C) scan (does not use the prepare module)."""
    rows = []
    for path in paths:
        lines = read_text(path).splitlines()
        index = 0
        while index < len(lines):
            count = int(lines[index].strip())
            comment = lines[index + 1]
            match = re.search(r"i\s*=\s*(\d+)[,\s]+time\s*=\s*([-0-9.]+)", comment)
            carbon = [float(value) for value in lines[index + 2 + 288].split()[1:4]]
            sulfur = [float(value) for value in lines[index + 2 + 292].split()[1:4]]
            delta = [sulfur[i] - carbon[i] for i in range(3)]
            delta = [delta[i] - box[i] * round(delta[i] / box[i]) for i in range(3)]
            rows.append((int(match.group(1)), float(match.group(2)),
                         (sum(value * value for value in delta)) ** 0.5))
            index += count + 2
    return rows


def check_source(manifest):
    source = manifest["source"]
    ok("seg2 source SHA256 matches the manifest",
       sha256(prep.SOURCE) == source["sha256"], sha256(prep.SOURCE))
    ok("seg2 source is the 1000-frame continuation (steps 10010-20000)",
       source["trajectory"].endswith("interface_r1_seg2/interface_r1-pos-1.xyz"))
    ok("first-segment trajectory is recorded with its SHA",
       manifest["twenty_ps_event_scan"]["trajectories"][0]["path"]
       .endswith("interface_r1_fix1/interface_r1-pos-1.xyz")
       and manifest["twenty_ps_event_scan"]["trajectories"][0]["sha256"]
       == sha256(prep.SOURCE_SEG1))
    ok("the scan covers both segments (1001 + 1000 = 2001 frames)",
       [item["frames"] for item in manifest["twenty_ps_event_scan"]["trajectories"]]
       == [1001, 1000]
       and manifest["twenty_ps_event_scan"]["frames_scanned"] == 2001)

    frame = bulk.find_event_frame_by_step(prep.SOURCE, prep.EVENT_STEP,
                                          prep.EVENT_TIME_FS, prep.TIME_TOLERANCE_FS)
    ok("exactly one frame matches step 10410 and time 10409.999160 +/- 0.01 fs",
       frame["step"] == 10410 and abs(frame["time_fs_in_file"] - 10410.0) < 1e-9)
    ok("the source comment prints three decimals, i.e. a different precision form",
       "time =    10410.000" in frame["comment"]
       and abs(frame["time_fs_in_file"] - prep.EVENT_TIME_FS) > 1e-4)
    ok("the numeric tolerance still selects that frame (no string comparison)",
       abs(frame["time_fs_in_file"] - prep.EVENT_TIME_FS) <= prep.TIME_TOLERANCE_FS)

    # Gates: a neighbouring step must not match, and zero matches must stop.
    for label, step, time_fs in (("previous step", 10400, 10409.999160),
                                 ("next step", 10420, 10409.999160),
                                 ("wrong time", 10410, 10419.999160)):
        try:
            bulk.find_event_frame_by_step(prep.SOURCE, step, time_fs,
                                          prep.TIME_TOLERANCE_FS)
            raise AssertionError("accepted a wrong (step, time) pair: %s" % label)
        except SystemExit:
            pass
    ok("a wrong step or a wrong time is rejected", True)

    elements, coords = bulk.parse_atoms(frame["atom_lines"])
    geometry = bulk.verify_structure(elements, coords, prep.BOX)
    ok("frame is 296 atoms H194 O100 C1 S1 in the fixed order",
       len(elements) == 296 and elements == prep.EXPECTED_ORDER and dict(
           (element, elements.count(element)) for element in ("H", "O", "C", "S"))
       == {"H": 194, "O": 100, "C": 1, "S": 1})
    ok("S is atom 293 and C is atom 289 (1-based)",
       elements.index("S") == 292 and elements.index("C") == 288)
    ok("S-C in the 15x15x30 cell is 1.808134 A",
       abs(geometry["sc_distance_angstrom"] - 1.808134) < 5e-4,
       "%.6f" % geometry["sc_distance_angstrom"])
    ok("carbonyl C=O %.4f A and C-H %.4f/%.4f A"
       % (geometry["co_distance_angstrom"], geometry["ch_distances_angstrom"][0],
          geometry["ch_distances_angstrom"][1]),
       abs(geometry["co_distance_angstrom"] - 1.342) < 5e-3
       and abs(geometry["ch_distances_angstrom"][0] - 1.0835) < 5e-3)
    ok("three S-O bonds are %.3f/%.3f/%.3f A (SO3 intact)"
       % tuple(geometry["so_distances_angstrom"]),
       all(abs(value - 1.52) < 0.08 for value in geometry["so_distances_angstrom"]))
    ok("no duplicate atoms and no severe overlap in the interface cell",
       geometry["min_atom_pair_distance_including_bonds_angstrom"] > 0.9
       and geometry["min_intermolecular_contact_angstrom"] > 1.2)
    ok("minimum intermolecular contact is %.4f A"
       % geometry["min_intermolecular_contact_angstrom"],
       abs(geometry["min_intermolecular_contact_angstrom"] - 1.4172) < 2e-3)

    # The same gate must fail with the bulk cube: the cells are not interchangeable.
    try:
        bulk.verify_structure(elements, coords, (14.2296, 14.2296, 14.2296))
        raise AssertionError("the bulk 14.2296 A cube was accepted for the interface frame")
    except SystemExit:
        ok("the bulk 14.2296 A cube is rejected for this interface frame", True)

    xyz_path = HERE / source["xyz"]["path"]
    xyz_lines = read_text(xyz_path).splitlines()
    ok("event XYZ is 298 lines with 296 atoms",
       len(xyz_lines) == 298 and int(xyz_lines[0]) == 296)
    ok("event XYZ atom lines are byte-identical to the trajectory",
       xyz_lines[2:] == frame["atom_lines"])
    ok("event XYZ SHA256 matches the manifest",
       sha256(xyz_path) == source["xyz"]["sha256"])
    ok("event XYZ comment is parseable and records both time forms",
       "i=10410" in xyz_lines[1] and "time=10409.999160" in xyz_lines[1]
       and "10410.000" in xyz_lines[1] and "verbatim" in xyz_lines[1])
    ok("manifest records the source line range, step and both times",
       source["frame"]["line_range_1based"][0] == frame["count_index"] + 1
       and source["frame"]["step"] == 10410
       and source["frame"]["time_fs"] == 10409.999160
       and source["frame"]["time_fs_printed_in_trajectory"] == 10410.0)

    scan = manifest["twenty_ps_event_scan"]
    rows = independent_sc_scan([prep.SOURCE_SEG1, prep.SOURCE], prep.BOX)
    short = [row for row in rows if row[2] < 2.0]
    ok("independent rescan also finds exactly six frames below 2.0 A",
       len(short) == 6 and scan["frames_below_threshold"] == 6, "%.0f" % len(short))
    ok("the six short contacts sit in one contiguous 10370-10420 window",
       [row[0] for row in short] == [10370, 10380, 10390, 10400, 10410, 10420]
       and scan["window_steps"] == [10370, 10420]
       and scan["single_window"] is True and scan["window_is_contiguous"] is True)
    ok("the window spans 10370-10420 fs with 10 fs spacing",
       scan["window_time_fs"] == [10370.0, 10420.0]
       and scan["window_step_spacing_fs"] == 10)
    shortest = min(rows, key=lambda row: row[2])
    ok("the 20 ps minimum is step 10410 at 1.808134 A",
       shortest[0] == 10410 and abs(shortest[2] - 1.808134) < 5e-4
       and scan["minimum"]["step"] == 10410)
    ok("manifest rows carry every short contact with step and time",
       [(row["step"], row["time_fs"]) for row in scan["below_threshold_rows"]]
       == [(row[0], row[1]) for row in short])
    ok("scan records the 2.0 A threshold and the interface cell",
       scan["threshold_angstrom"] == 2.0 and scan["box_angstrom"] == [15.0, 15.0, 30.0])

    colvar = manifest["colvar_cross_check"]
    ok("COLVAR cross-check finds exactly one row at 10409.999160 fs",
       colvar["matches"] == 1 and abs(colvar["time_fs"] - 10409.999160) < 1e-9)
    ok("COLVAR d_sc agrees with the extracted frame to 5e-4 A",
       abs(colvar["d_sc_angstrom"] - geometry["sc_distance_angstrom"]) < 5e-4)
    ok("manifest declares the interface cell and its whitelist class",
       manifest["box_angstrom"] == [15.0, 15.0, 30.0]
       and manifest["box_class"] == "interface")
    ok("manifest records no verdict", manifest.get("verdict") is None)
    return geometry


def check_manifest(manifest):
    ok("manifest schema", manifest.get("schema")
       == "li2025_direct_interface_event_unbiased/v1")
    cases = manifest["cases"]
    ok("three replicas are defined",
       [case["name"] for case in cases] == ["interface_event_unbiased_s1",
                                            "interface_event_unbiased_s2",
                                            "interface_event_unbiased_s3"])
    ok("seeds are 20261112/13/14",
       [case["seed"] for case in cases] == [20261112, 20261113, 20261114])
    ok("every replica is unbiased, interface-class and has no FES",
       all(case["unbiased"] is True and case["box_class"] == "interface"
           and case["fes"]["status"] == "not_applicable"
           and case["box_angstrom"] == [15.0, 15.0, 30.0] for case in cases))
    ok("each replica runs 3000 steps / 301 frames",
       all(case["steps"] == 3000 and case["expected_frames"] == 301 for case in cases))
    ok("outputs follow the interface_r1_unbiased project",
       all(case["outputs"]["trajectory"] == "interface_r1_unbiased-pos-1.xyz"
           and case["outputs"]["energy"] == "interface_r1_unbiased-1.ener"
           for case in cases))
    ok("velocity strategy is explicitly fresh, not inherited",
       manifest["velocity_strategy"]["reused_metad_velocities"] is False
       and "no &EXT_RESTART" in manifest["velocity_strategy"]["how"])
    scope = " ".join(manifest["scope"]["what_this_is_not"]).lower()
    ok("scope disclaims spontaneous formation, stability and free energy",
       "spontaneous" in scope and "stable" in scope and "free energy" in scope)
    ok("production records the interface cell and q = -2 with no counterion",
       manifest["production"]["cell"] == "ABC 15.0 15.0 30.0"
       and manifest["production"]["charge"] == -2
       and manifest["production"]["counterion"] is None)
    return cases


def check_inputs(cases, geometry):
    reference = read_text(prep.REFERENCE_INPUT)
    seeds = [case["seed"] for case in cases]
    for case in cases:
        text = read_text(HERE / "inputs" / (case["name"] + ".inp"))
        upper = text.upper()
        for forbidden in ("PLUMED", "METADYN", "FREE_ENERGY", "USE_PLUMED",
                          "PLUMED_INPUT_FILE", "EXT_RESTART", "VELOCITY"):
            ok("%s: input has no %s" % (case["name"], forbidden), forbidden not in upper)
        ok("%s: STEPS is 3000" % case["name"],
           re.search(r"^\s*STEPS 3000\s*$", text, re.MULTILINE) is not None)
        ok("%s: SEED is inside &GLOBAL" % case["name"],
           text.index("  SEED %d\n" % case["seed"]) < text.index("&END GLOBAL"))
        ok("%s: cell is ABC 15.0 15.0 30.0 and no 14.2296 appears" % case["name"],
           "ABC 15.0 15.0 30.0" in text and "14.2296" not in text)
        for token in ("ENSEMBLE NVT", "TYPE CSVR", "TIMECON 100", "TIMESTEP 1.0",
                      "TEMPERATURE 298.15", "&XC_FUNCTIONAL BLYP", "TYPE DFTD3(BJ)",
                      "REFERENCE_FUNCTIONAL BLYP", "CHARGE -2", "MULTIPLICITY 1",
                      "CUTOFF 300", "REL_CUTOFF 40", "DZVP-MOLOPT-SR-GTH",
                      "GTH-BLYP-q1", "GTH-BLYP-q6", "GTH-BLYP-q4",
                      "COORD_FILE_NAME interface_r1_event10410.xyz",
                      "MD 10", "MD 100", "EPS_SCF 1.0E-6", "MINIMIZER DIIS"):
            ok("%s: keeps %r" % (case["name"], token), token in text)

        removed, _block = bulk.strip_free_energy(reference)
        diff_lines = []
        for line in difflib.unified_diff(removed.splitlines(), text.splitlines(),
                                         lineterm="", n=0):
            if line.startswith(("+++", "---")):
                continue
            if line.startswith(("+", "-")):
                diff_lines.append(line)
        expected_removed = {
            "  PROJECT interface_r1",
            "    STEPS 10000",
            "      COORD_FILE_NAME interface_r1.xyz",
        }
        expected_added = {
            "  PROJECT interface_r1_unbiased",
            "    STEPS 3000",
            "      COORD_FILE_NAME interface_r1_event10410.xyz",
            "  SEED %d" % case["seed"],
        }
        actual_removed = set(line[1:] for line in diff_lines if line.startswith("-"))
        actual_added = set(line[1:] for line in diff_lines if line.startswith("+"))
        ok("%s: only the intended lines differ from interface_r1.inp" % case["name"],
           actual_removed == expected_removed and actual_added == expected_added,
           "removed=%s added=%s" % (sorted(actual_removed), sorted(actual_added)))

    texts = [read_text(HERE / "inputs" / (case["name"] + ".inp")) for case in cases]
    normalised = [re.sub(r"^  SEED \d+\n", "  SEED <seed>\n", text, flags=re.MULTILINE)
                  for text in texts]
    ok("the three replica inputs differ only in the SEED line",
       len(set(normalised)) == 1 and len(set(seeds)) == 3)
    ok("each input stages the shared interface event XYZ",
       all(any(entry["source"] == "interface_r1_event10410.xyz" for entry in case["stage"])
           for case in cases))
    ok("allocated stage SHAs match the files on disk",
       all(sha256(HERE / "inputs" / entry["source"]) == entry["sha256"]
           for case in cases for entry in case["stage"]))
    ok("inputs record the frame S-C in the manifest case block",
       all(abs(case["sc_distance_angstrom"] - geometry["sc_distance_angstrom"]) < 1e-9
           for case in cases))


def check_build_gates():
    reference = read_text(prep.REFERENCE_INPUT)
    try:
        bulk.build_unbiased_input(reference, 1, prep.FRAME_STEM,
                                  project=prep.PROJECT, steps=prep.STEPS,
                                  reference_project="interface_r1",
                                  reference_coord="interface_r1.xyz",
                                  expected_cell="ABC 15.0 15.0 30.0")
        ok("the interface reference builds when the expected cell is stated", True)
    except SystemExit:
        ok("the interface reference builds when the expected cell is stated", False,
           "build refused")
    try:
        bulk.build_unbiased_input(reference, 1, prep.FRAME_STEM,
                                  project=prep.PROJECT, steps=prep.STEPS,
                                  reference_project="interface_r1",
                                  reference_coord="interface_r1.xyz",
                                  expected_cell="ABC 14.2296 14.2296 14.2296")
        raise AssertionError("accepted a bulk cell for the interface input")
    except SystemExit:
        ok("a build that asks for the bulk cell is refused", True)


def check_wrappers():
    compute = [HERE / ("submit_interface_event_unbiased_s%d.slurm" % index)
               for index in (1, 2, 3)]
    gpu = [HERE / ("submit_interface_event_unbiased_s%d_gpu.slurm" % index)
           for index in (1, 2, 3)]
    workdirs = set()
    job_names = set()
    logs = set()
    for path in compute + gpu:
        is_gpu = path.name.endswith("_gpu.slurm")
        ok("%s exists and passes bash -n" % path.name,
           path.is_file() and subprocess.call(["bash", "-n", str(path)]) == 0)
        text = read_text(path)
        ok("%s: 64 ranks, 32G, 48 h, mail END,FAIL" % path.name,
           "#SBATCH --ntasks=64" in text and "#SBATCH --mem=32G" in text
           and "#SBATCH --time=48:00:00" in text
           and "#SBATCH --mail-type=END,FAIL" in text)
        sbatch_lines = [line for line in text.splitlines()
                        if line.startswith("#SBATCH")]
        ok("%s: no --gres and one thread per rank" % path.name,
           all("--gres" not in line for line in sbatch_lines)
           and "OMP_NUM_THREADS=1" in text)
        ok("%s: partition is %s" % (path.name, "gpu" if is_gpu else "compute"),
           ("#SBATCH --partition=gpu" in text) if is_gpu
           else ("#SBATCH --partition=compute" in text))
        if is_gpu:
            ok("%s: requests one CPU per rank, no GPU" % path.name,
               "#SBATCH --cpus-per-task=1" in text
               and all("gpu:" not in line for line in sbatch_lines))
        else:
            ok("%s: requests no --cpus-per-task" % path.name,
               "--cpus-per-task" not in text)
        job = re.search(r"^#SBATCH --job-name=(\S+)$", text, re.MULTILINE)
        work = re.search(r"^WORK=\$BASE/(\S+)$", text, re.MULTILINE)
        ok("%s: distinct job name and workdir" % path.name,
           job is not None and work is not None)
        if job and work:
            job_names.add(job.group(1))
            workdirs.add(work.group(1))
        logs.add(re.search(r"^#SBATCH --output=(\S+)$", text, re.MULTILINE).group(1))
        ok("%s: uses the interface manifest and the matching case" % path.name,
           "--manifest %s" % MANIFEST_NAME in text
           and ("--case " + work.group(1)) in text)
        executable = "\n".join(line for line in text.splitlines()
                               if not line.lstrip().startswith("#"))
        ok("%s: analyzes with the interface cell, never the bulk one" % path.name,
           "analyze_event_unbiased.py" in executable
           and "--box 15.0,15.0,30.0" in executable
           and "14.2296" not in executable)
        ok("%s: reads the interface_r1_unbiased outputs" % path.name,
           "interface_r1_unbiased-pos-1.xyz" in text
           and "interface_r1_unbiased-1.ener" in text)
    ok("six distinct workdirs",
       workdirs == {"interface_event_unbiased_s1", "interface_event_unbiased_s2",
                    "interface_event_unbiased_s3"})
    ok("six distinct job names", len(job_names) == 6, str(sorted(job_names)))
    ok("six distinct log prefixes", len(logs) == 6)


def check_runner(cases):
    scratch = tempfile.mkdtemp(prefix="li2025_iface_unbiased_selftest_")
    try:
        for case in cases:
            name = case["name"]
            workdir = Path(scratch) / name
            code, output = run_runner([
                "--manifest", MANIFEST_NAME, "--case", name,
                "--workdir", str(workdir), "--workdir-root", scratch,
                "--ntasks", "64", "--dry-run"])
            ok("%s: dry-run succeeds" % name, code == 0, output.strip()[-200:])
            payload = json.loads(output)
            ok("%s: dry-run stages two files and verifies hashes" % name,
               payload["dry_run"] and len(payload["staged"]) == 2)
            ok("%s: dry-run reports the interface cell, not the bulk cube" % name,
               payload["box_class"] == "interface"
               and payload["box_angstrom"] == [15.0, 15.0, 30.0]
               and payload["input_cell_angstrom"] == [15.0, 15.0, 30.0])
            ok("%s: dry-run reports the seed and the 3000-step budget" % name,
               payload["unbiased"]["seed"] == case["seed"]
               and payload["unbiased"]["steps"] == 3000
               and payload["unbiased"]["expected_frames"] == 301)
            ok("%s: staged input equals the generated input" % name,
               sha256(workdir / (name + ".inp"))
               == sha256(HERE / "inputs" / (name + ".inp")))
            ok("%s: no PLUMED file is staged" % name,
               not (workdir / (name + ".plumed")).exists())

        bulk_case = cases[0]

        # Hash mismatch.
        workdir = Path(scratch) / "tamper"
        code, _ = run_runner([
            "--manifest", MANIFEST_NAME, "--case", bulk_case["name"],
            "--workdir", str(workdir), "--workdir-root", scratch, "--dry-run"])
        ok("first dry-run into the tamper workdir succeeds", code == 0)
        xyz = workdir / "interface_r1_event10410.xyz"
        xyz.write_text(xyz.read_text(encoding="utf-8") + "\n", encoding="utf-8")
        code, _ = run_runner([
            "--manifest", MANIFEST_NAME, "--case", bulk_case["name"],
            "--workdir", str(workdir), "--workdir-root", scratch, "--dry-run"])
        ok("runner refuses a tampered event XYZ (exit 4)", code == 4)

        # Stale output.
        stale = Path(scratch) / "stale"
        stale.mkdir()
        (stale / (bulk_case["name"] + ".out")).write_text("old\n", encoding="utf-8")
        code, _ = run_runner([
            "--manifest", MANIFEST_NAME, "--case", bulk_case["name"],
            "--workdir", str(stale), "--workdir-root", scratch, "--dry-run"])
        ok("runner refuses a stale output workdir (exit 7)", code == 7)

        # Unsafe workdirs.
        code, _ = run_runner([
            "--manifest", MANIFEST_NAME, "--case", bulk_case["name"],
            "--workdir", str(Path(scratch) / ".." / "escape"),
            "--workdir-root", scratch, "--dry-run"])
        ok("runner refuses '..' in the workdir (exit 5)", code == 5)
        code, _ = run_runner([
            "--manifest", MANIFEST_NAME, "--case", bulk_case["name"],
            "--workdir", str(Path(scratch) / "nested" / "deep"),
            "--workdir-root", scratch, "--dry-run"])
        ok("runner refuses a nested workdir (exit 5)", code == 5)

        # Case/manifest whitelist.
        code, output = run_runner([
            "--manifest", MANIFEST_NAME, "--case", "event_unbiased_s1",
            "--workdir", str(Path(scratch) / "wrongcase"),
            "--workdir-root", scratch, "--dry-run"])
        ok("runner refuses a bulk case name for the interface manifest",
           code != 0 and "not registered" in output)
        bogus_dir = Path(scratch) / "bogus_manifest"
        bogus_dir.mkdir()
        bogus = bogus_dir / "unregistered_manifest.json"
        bogus.write_text(read_text(MANIFEST), encoding="utf-8")
        code, output = run_runner([
            "--manifest", str(bogus), "--case", bulk_case["name"],
            "--workdir", str(Path(scratch) / "bogus"),
            "--workdir-root", scratch, "--dry-run"])
        ok("runner refuses an unregistered manifest file name",
           code != 0 and "unregistered manifest" in output)

        # box_class vs box_angstrom (manifest side).
        tampered = json.loads(read_text(MANIFEST))
        for entry in tampered["cases"]:
            entry["box_angstrom"] = [14.2296, 14.2296, 14.2296]
        tampered_dir = Path(scratch) / "tampered_manifest"
        tampered_dir.mkdir()
        tampered_path = tampered_dir / MANIFEST_NAME
        tampered_path.write_text(json.dumps(tampered, indent=2), encoding="utf-8")
        code, output = run_runner([
            "--manifest", str(tampered_path), "--case", bulk_case["name"],
            "--workdir", str(Path(scratch) / "badbox"),
            "--workdir-root", scratch, "--dry-run"])
        ok("runner refuses an interface case whose box_angstrom is the bulk cube (exit 4)",
           code == 4 and "box_class" in output)

        # box_class vs the staged ABC line (input side).
        cell_dir = Path(scratch) / "cell_gate"
        cell_dir.mkdir()
        (cell_dir / "fake.xyz").write_text("1\nx\nH 0 0 0\n", encoding="utf-8")
        (cell_dir / "fake.inp").write_text(
            "&SUBSYS\n  &CELL\n    ABC 14.2296 14.2296 14.2296\n  &END CELL\n"
            "  &TOPOLOGY\n    COORD_FILE_NAME fake.xyz\n  &END TOPOLOGY\n&END SUBSYS\n",
            encoding="utf-8")
        entries = [{"target": "fake.inp"}, {"target": "fake.xyz"}]
        try:
            runner.check_case_inputs(cell_dir, "fake", {"box_class": "interface"}, entries)
            raise AssertionError("accepted the bulk cell for an interface case")
        except RuntimeError as error:
            ok("runner refuses a staged bulk ABC line for an interface case",
               "whitelist cell" in str(error))

        (cell_dir / "fake.inp").write_text(
            "&SUBSYS\n  &CELL\n    ABC 15.0 15.0 30.0\n  &END CELL\n"
            "  &TOPOLOGY\n    COORD_FILE_NAME other.xyz\n  &END TOPOLOGY\n&END SUBSYS\n",
            encoding="utf-8")
        try:
            runner.check_case_inputs(cell_dir, "fake", {"box_class": "interface"}, entries)
            raise AssertionError("accepted an unstaged COORD_FILE_NAME")
        except RuntimeError as error:
            ok("runner refuses a COORD_FILE_NAME the staging plan does not provide",
               "staging plan does not provide" in str(error))

        (cell_dir / "fake.inp").write_text(
            "&SUBSYS\n  &CELL\n    ABC 15.0 15.0 30.0\n  &END CELL\n"
            "  &TOPOLOGY\n    COORD_FILE_NAME fake.xyz\n  &END TOPOLOGY\n&END SUBSYS\n",
            encoding="utf-8")
        ok("runner accepts the interface cell and a staged coordinate file",
           runner.check_case_inputs(cell_dir, "fake", {"box_class": "interface"},
                                    entries) == [15.0, 15.0, 30.0])

        # End-to-end with fake CP2K output: complete (301 frames) and truncated.
        for label, last_step, frames, expected_code in (
                ("complete", 3000, 301, 0), ("short", 2999, 301, 6)):
            work = Path(scratch) / ("fake_" + label)

            def fake_cp2k(command, cwd=None, env=None, stdout=None, stderr=None,
                          last_step=last_step, frames=frames):
                Path(cwd, "interface_event_unbiased_s1.out").write_text("\n".join([
                    " MD| Step number 1",
                    " SCF run converged in 8 steps",
                    " ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]:  -1742.10",
                    " MD| Step number %d" % last_step,
                    " SCF run converged in 9 steps",
                    " ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]:  -1742.11",
                    " PROGRAM ENDED AT   test",
                ]) + "\n", encoding="utf-8")
                body = []
                for index in range(frames):
                    step = index * 10
                    body.append("3\n i = %d, time = %d.000, E = -1.0\n"
                                "H 0 0 0\nH 0 0 1\nO 0 0 2\n" % (step, step))
                Path(cwd, "interface_r1_unbiased-pos-1.xyz").write_text(
                    "".join(body), encoding="utf-8")
                Path(cwd, "interface_r1_unbiased-1.ener").write_text(
                    "#     Step Nr.          Time[fs]        Kin.[a.u.]        "
                    "Temp[K]            Pot.[a.u.]        Cons Qty[a.u.]        UsedTime[s]\n"
                    "         0            0.000000         0.4       298.150000     "
                    "-1742.1       -1741.7       0.0\n"
                    "      %d          %d.000000         0.4       298.150000     "
                    "-1742.1       -1741.7       0.0\n" % (last_step, last_step),
                    encoding="utf-8")
                return 0

            original_call = runner.subprocess.call
            try:
                runner.subprocess.call = fake_cp2k
                code = runner.run("interface_event_unbiased_s1", str(work), 64, False,
                                  MANIFEST_NAME, None, scratch)
            finally:
                runner.subprocess.call = original_call
            record = json.loads(read_text(work / "run_record.json"))
            if expected_code == 0:
                ok("end-to-end complete run exits 0 and records COMPLETED",
                   code == 0 and record["status"] == "COMPLETED"
                   and record["unbiased"]["step_check"] == "PASS"
                   and record["unbiased"]["trajectory_check"] == "PASS"
                   and record["unbiased"]["trajectory_frames"] == 301
                   and record["fes"]["status"] == "not_applicable"
                   and "metad" not in record)
                ok("run record keeps the interface cell and box class",
                   record["box_class"] == "interface"
                   and record["box_angstrom"] == [15.0, 15.0, 30.0]
                   and record["input_cell_angstrom"] == [15.0, 15.0, 30.0])
            else:
                ok("end-to-end truncated run exits 6 and records FAILED_STEP_RANGE",
                   code == 6 and record["status"] == "FAILED_STEP_RANGE"
                   and record["unbiased"]["step_check"] == "FAIL")
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def check_analyzer():
    scratch = tempfile.mkdtemp(prefix="li2025_iface_analyzer_")
    try:
        xyz = HERE / "inputs" / "interface_r1_event10410.xyz"
        out = Path(scratch) / "single.json"
        code, output = run_process([
            sys.executable, str(HERE / "analyze_event_unbiased.py"),
            "--trajectory", str(xyz), "--label", "interface_event",
            "--box", INTERFACE_BOX, "--out", str(out)])
        ok("analyzer runs on the extracted interface frame", code == 0, output.strip()[-200:])
        data = json.loads(read_text(out))
        ok("analyzer reports one frame with S-C 1.808134 A",
           data["frames"]["count"] == 1
           and abs(data["per_frame"][0]["sc"] - 1.808134) < 5e-4)
        ok("analyzer records the interface cell it was given",
           data["box_angstrom"] == [15.0, 15.0, 30.0])
        ok("analyzer asserts no verdict",
           data["interpretation"]["verdict"] is None
           and "unbiased" in data["interpretation"]["what_this_is"].lower())
        ok("analyzer disclaims free energy and spontaneous formation",
           any("free energy" in item for item in data["interpretation"]["what_this_is_not"])
           and any("spontaneous" in item
                   for item in data["interpretation"]["what_this_is_not"]))
        tsv = out.with_name(out.stem + "_FRAMES.tsv")
        ok("analyzer writes a per-frame TSV with header + one row",
           tsv.is_file() and len(read_text(tsv).splitlines()) == 2)

        code, _ = run_process([
            sys.executable, str(HERE / "analyze_event_unbiased.py"),
            "--trajectory", str(xyz), "--out", str(Path(scratch) / "nobox.json")])
        ok("analyzer refuses to run without an explicit --box", code != 0)

        code, _ = run_process([
            sys.executable, str(HERE / "analyze_event_unbiased.py"),
            "--trajectory", str(xyz), "--box", "15.0,15.0",
            "--out", str(Path(scratch) / "shortbox.json")])
        ok("analyzer rejects a two-component --box", code != 0)

        # A synthetic 3-frame interface-box trajectory with a clear dissociation.
        frames = []
        for index, sc in enumerate((1.8, 2.6, 3.4)):
            body = []
            for water in range(96):
                x = (water % 12) * 1.1 + 0.5
                y = 2.0 + (water // 12) * 0.4
                body += [("H", x, y, 2.0), ("H", x + 0.98, y, 2.0),
                         ("O", x + 0.49, y + 0.85, 2.0)]
            cx, cy, cz = 2.0, 12.0, 8.0
            body += [("C", cx, cy, cz), ("O", cx + 1.2, cy, cz),
                     ("H", cx - 1.1, cy, cz), ("H", cx, cy + 1.1, cz)]
            sx = cx + sc
            body += [("S", sx, cy, cz), ("O", sx + 1.5, cy, cz),
                     ("O", sx, cy + 1.5, cz), ("O", sx, cy, cz + 1.5)]
            frames.append(body)
        trajectory = Path(scratch) / "interface_synthetic.xyz"
        lines = []
        for index, body in enumerate(frames):
            lines.append(str(len(body)))
            lines.append(" i = %d, time = %d.000, E = -1.0" % (index, index))
            for element, x, y, z in body:
                lines.append("%-2s %14.9f %14.9f %14.9f" % (element, x, y, z))
        trajectory.write_text("\n".join(lines) + "\n", encoding="utf-8")
        synth_out = Path(scratch) / "interface_synthetic.json"
        code, _ = run_process([
            sys.executable, str(HERE / "analyze_event_unbiased.py"),
            "--trajectory", str(trajectory), "--box", INTERFACE_BOX,
            "--out", str(synth_out)])
        synth = json.loads(read_text(synth_out)) if code == 0 else None
        ok("analyzer tracks S-C 1.8 -> 3.4 A in the interface cell",
           code == 0 and abs(synth["summary"]["sc"]["first"] - 1.8) < 1e-6
           and abs(synth["summary"]["sc"]["last"] - 3.4) < 1e-6
           and synth["box_angstrom"] == [15.0, 15.0, 30.0])
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def check_python36_and_shell():
    py_files = sorted(HERE.glob("*.py"))
    ok("the interface event scripts are present",
       {"prepare_interface_event_unbiased.py", "analyze_event_unbiased.py",
        "selftest_interface_event_unbiased.py"}.issubset(
            set(path.name for path in py_files)))
    for path in py_files:
        source = read_text(path)
        try:
            ast.parse(source, filename=str(path), feature_version=(3, 6))
            compile(source, str(path), "exec")
        except SyntaxError as error:
            raise AssertionError("%s is not Python 3.6 syntax: %s" % (path.name, error))
        capture_keyword = "capture" + "_output"
        text_keyword = "te" + "xt"
        for token in (capture_keyword, "missing" + "_ok",
                      "remove" + "prefix", "remove" + "suffix"):
            ok("%s has no %s" % (path.name, token), token not in source)
        pattern = (r"subprocess\.(run|Popen|call|check_call|check_output)"
                   r"\([^)]*?\b(" + text_keyword + "|" + capture_keyword + r")\s*=")
        ok("%s has no subprocess text/capture keywords" % path.name,
           re.search(pattern, source, re.DOTALL) is None)
        ok("%s has no walrus operator" % path.name, (":" + "=") not in source)
    for path in sorted(HERE.glob("*.slurm")):
        ok("%s passes bash -n" % path.name,
           subprocess.call(["bash", "-n", str(path)]) == 0)


def check_scope(manifest):
    readme = read_text(HERE / "README.md").lower()
    for token in ("interface_event_unbiased", "15 × 15 × 30", "10410",
                  "not a free energy", "gpu"):
        ok("README documents %r" % token, token in readme)
    leftovers = [path.name for path in HERE.rglob("*")
                 if path.suffix.lower() in (".hills", ".colvar")]
    ok("no HILLS/COLVAR files in the campaign directory", not leftovers, str(leftovers))
    ok("the manifest states no FES was produced",
       all(case["fes"]["status"] == "not_applicable" for case in manifest["cases"]))


def main():
    manifest = json.loads(read_text(MANIFEST))
    geometry = check_source(manifest)
    cases = check_manifest(manifest)
    check_inputs(cases, geometry)
    check_build_gates()
    check_wrappers()
    check_runner(cases)
    check_analyzer()
    check_python36_and_shell()
    check_scope(manifest)
    print("SELFTEST_INTERFACE_EVENT_UNBIASED_OK: %d checks passed (source, manifest, "
          "inputs, build gates, wrappers, runner, analyzer, py36, scope)" % CHECKS["count"])
    return 0


if __name__ == "__main__":
    sys.exit(main())

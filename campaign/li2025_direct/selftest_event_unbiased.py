#!/usr/bin/env python3
"""Self-test for the unbiased Li 2025 event replicas.

Uses the real pulled event trajectory for source/integration checks and
synthetic trajectories for the analyzer and the hard gates. No cluster, no CP2K.

Covered:
  * source frame uniqueness at time = 12270.000 fs and verbatim coordinate copy,
  * composition, atom order, S-C/C=O/C-H/S-O, contacts and event provenance,
  * unbiased input: no PLUMED/METADYN/FREE_ENERGY/EXT_RESTART, STEPS 3000,
    fresh-velocity policy (no VELOCITY, no EXT_RESTART), and three replicas that
    differ only in &GLOBAL SEED,
  * three wrappers: 64 ranks, 32G, 48 h, mail END/FAIL, no --gres, distinct
    workdirs/logs, analyzer post-step,
  * runner: dry-run staging, hash mismatch (exit 4), stale output (exit 7),
    unsafe workdir (exit 5), complete and truncated fake runs (exit 0 / 6),
  * analyzer: per-frame metrics, <2 A run, first sustained dissociation,
    temperature coverage, no verdict, TSV output, bad atom order rejected,
  * Python 3.6 syntax/API and `bash -n` for every script in the directory.

Run: python3 selftest_event_unbiased.py
Exit 0 and a final SELFTEST_OK line on success, 1 on any failure.
"""
import sys

sys.dont_write_bytecode = True

import ast
import hashlib
import json
import math
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))

import prepare_event_unbiased as prep
import analyze_event_unbiased as analyzer
import run_cp2k as runner

MANIFEST = HERE / "event_unbiased_manifest.json"
BULK_BOX = "14.2296,14.2296,14.2296"
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


def check_source(manifest):
    source = manifest["source"]
    ok("source SHA256 is the recorded value",
       sha256(prep.SOURCE) == source["sha256"], sha256(prep.SOURCE))

    frame = prep.find_event_frame(prep.SOURCE, prep.EVENT_TIME_FS)
    ok("exactly one frame matches time 12270.000 fs", frame["step"] == 12270)

    with open(str(prep.SOURCE), "r") as handle:
        lines = handle.read().splitlines()
    matches = [line for line in lines
               if line.lstrip().startswith("i =") and "time =    12270.000" in line]
    ok("manual rescan also finds exactly one 12270.000 fs frame", len(matches) == 1)

    elements, coords = prep.parse_atoms(frame["atom_lines"])
    geometry = prep.verify_structure(elements, coords)
    ok("frame is 296 atoms H194 O100 C1 S1 in the fixed order",
       len(elements) == 296 and elements == prep.EXPECTED_ORDER)
    ok("S is atom 293 and C is atom 289 (1-based)",
       elements.index("S") == 292 and elements.index("C") == 288)
    ok("S-C is 1.778073 A",
       abs(geometry["sc_distance_angstrom"] - 1.778073) < 5.0e-4,
       "%.6f" % geometry["sc_distance_angstrom"])
    ok("carbonyl C=O %.4f A and C-H %.4f/%.4f A"
       % (geometry["co_distance_angstrom"], geometry["ch_distances_angstrom"][0],
          geometry["ch_distances_angstrom"][1]),
       abs(geometry["co_distance_angstrom"] - 1.387573) < 2.0e-3
       and abs(geometry["ch_distances_angstrom"][0] - 1.069856) < 2.0e-3)
    ok("three S-O bonds are %.3f/%.3f/%.3f A"
       % tuple(geometry["so_distances_angstrom"]),
       all(abs(value - 1.5) < 0.06 for value in geometry["so_distances_angstrom"]))
    ok("no duplicate atoms and no severe overlap",
       geometry["min_atom_pair_distance_including_bonds_angstrom"] > 0.9
       and geometry["min_intermolecular_contact_angstrom"] > 1.2)
    ok("minimum intermolecular contact is the O290-H268 hydrogen bond",
       abs(geometry["min_intermolecular_contact_angstrom"] - 1.478026) < 2.0e-3
       and {geometry["min_intermolecular_contact_atoms_1based"][0],
            geometry["min_intermolecular_contact_atoms_1based"][2]} == {268, 290})

    xyz_path = HERE / source["xyz"]["path"]
    xyz_lines = read_text(xyz_path).splitlines()
    ok("event XYZ is 298 lines with 296 atoms", len(xyz_lines) == 298
       and int(xyz_lines[0]) == 296)
    ok("event XYZ atom lines are byte-identical to the trajectory",
       xyz_lines[2:] == frame["atom_lines"])
    ok("event XYZ SHA256 matches the manifest",
       sha256(xyz_path) == source["xyz"]["sha256"])
    ok("event XYZ comment records frame provenance",
       "i=12270" in xyz_lines[1] and "1.778073" in xyz_lines[1]
       and "verbatim" in xyz_lines[1])
    ok("manifest records the source line range and step/time",
       source["frame"]["line_range_1based"][0] == frame["count_index"] + 1
       and source["frame"]["time_fs"] == 12270.0
       and source["frame"]["step"] == 12270)

    event = manifest["event"]
    ok("event analysis JSON is recorded with a SHA",
       event.get("path") is not None and event.get("sha256") == sha256(prep.EVENT_ANALYSIS))
    ok("event analysis S-C at 12270 agrees with the frame",
       abs(event["sc_at_event_from_analysis_angstrom"]
           - geometry["sc_distance_angstrom"]) < 1.0e-6)
    ok("12270 is the minimum inside the 12220-12320 fs event window",
       abs(event["event_window_min_sc_angstrom"] - 1.778073) < 5.0e-4)
    ok("manifest records that 12270 is not the trajectory-wide minimum",
       event["requested_frame_is_global_min"] is False
       and abs(event["trajectory_min_sc_angstrom"] - 1.7544) < 5.0e-4
       and event["trajectory_min_sc_step"] == 13320
       and "NOT used here" in event["note"])
    return geometry


def check_manifest(manifest):
    ok("manifest schema", manifest.get("schema") == "li2025_direct_event_unbiased/v1")
    cases = manifest["cases"]
    ok("three replicas are defined",
       [case["name"] for case in cases]
       == ["event_unbiased_s1", "event_unbiased_s2", "event_unbiased_s3"])
    ok("seeds are 20261109/10/11",
       [case["seed"] for case in cases] == [20261109, 20261110, 20261111])
    ok("all replicas are labelled unbiased with no FES",
       all(case["unbiased"] is True and case["fes"]["status"] == "not_applicable"
           for case in cases))
    ok("each replica runs 3000 steps / 301 frames",
       all(case["steps"] == 3000 and case["expected_frames"] == 301 for case in cases))
    ok("velocity strategy is explicitly fresh, not inherited",
       manifest["velocity_strategy"]["reused_metad_velocities"] is False
       and "no &EXT_RESTART" in manifest["velocity_strategy"]["how"])
    scope = " ".join(manifest["scope"]["what_this_is_not"]).lower()
    ok("scope disclaims spontaneous formation and free energy",
       "spontaneous" in scope and "free energy" in scope)
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
        for token in ("ENSEMBLE NVT", "TYPE CSVR", "TIMECON 100", "TIMESTEP 1.0",
                      "TEMPERATURE 298.15", "&XC_FUNCTIONAL BLYP", "TYPE DFTD3(BJ)",
                      "REFERENCE_FUNCTIONAL BLYP", "CHARGE -2", "MULTIPLICITY 1",
                      "ABC 14.2296 14.2296 14.2296", "CUTOFF 300", "REL_CUTOFF 40",
                      "DZVP-MOLOPT-SR-GTH", "GTH-BLYP-q1", "GTH-BLYP-q6",
                      "COORD_FILE_NAME bulk_r1_event12270.xyz",
                      "MD 10", "MD 100", "EPS_SCF 1.0E-6", "MINIMIZER DIIS"):
            ok("%s: keeps %r" % (case["name"], token), token in text)

        # Independent text-level check: only the intended edits.
        removed, _block = prep.strip_free_energy(reference)
        diff_lines = []
        import difflib
        for line in difflib.unified_diff(removed.splitlines(), text.splitlines(),
                                         lineterm="", n=0):
            if line.startswith(("+++", "---")):
                continue
            if line.startswith(("+", "-")):
                diff_lines.append(line)
        expected_removed = {
            "  PROJECT bulk_r1",
            "    STEPS 10000",
            "      COORD_FILE_NAME bulk_r1.xyz",
        }
        expected_added = {
            "  PROJECT bulk_r1_unbiased",
            "    STEPS 3000",
            "      COORD_FILE_NAME bulk_r1_event12270.xyz",
            "  SEED %d" % case["seed"],
        }
        actual_removed = set(line[1:] for line in diff_lines if line.startswith("-"))
        actual_added = set(line[1:] for line in diff_lines if line.startswith("+"))
        ok("%s: only the intended lines differ from the reference" % case["name"],
           actual_removed == expected_removed and actual_added == expected_added,
           "removed=%s added=%s" % (sorted(actual_removed), sorted(actual_added)))

    texts = [read_text(HERE / "inputs" / (case["name"] + ".inp")) for case in cases]
    normalised = [re.sub(r"^  SEED \d+\n", "  SEED <seed>\n", text, flags=re.MULTILINE)
                  for text in texts]
    ok("the three replica inputs differ only in the SEED line",
       len(set(normalised)) == 1 and len(set(seeds)) == 3)
    ok("each input stages the shared event XYZ",
       all(any(entry["source"] == "bulk_r1_event12270.xyz" for entry in case["stage"])
           for case in cases))
    ok("allocated stage SHAs match the files on disk",
       all(sha256(HERE / "inputs" / entry["source"]) == entry["sha256"]
           for case in cases for entry in case["stage"]))


def check_wrappers():
    wrappers = [HERE / ("submit_event_unbiased_s%d.slurm" % index) for index in (1, 2, 3)]
    workdirs = set()
    job_names = set()
    logs = set()
    for wrapper in wrappers:
        ok("%s exists and passes bash -n" % wrapper.name,
           wrapper.is_file() and subprocess.call(["bash", "-n", str(wrapper)]) == 0)
        text = read_text(wrapper)
        ok("%s: 64 ranks, 32G, 48 h, mail END,FAIL" % wrapper.name,
           "#SBATCH --ntasks=64" in text and "#SBATCH --mem=32G" in text
           and "#SBATCH --time=48:00:00" in text
           and "#SBATCH --mail-type=END,FAIL" in text)
        ok("%s: no GPU request and one thread per rank" % wrapper.name,
           "--gres" not in text and "--cpus-per-task" not in text
           and "OMP_NUM_THREADS=1" in text)
        job = re.search(r"^#SBATCH --job-name=(\S+)$", text, re.MULTILINE)
        work = re.search(r"^WORK=\$BASE/(\S+)$", text, re.MULTILINE)
        ok("%s: distinct job name and workdir" % wrapper.name,
           job is not None and work is not None)
        if job and work:
            job_names.add(job.group(1))
            workdirs.add(work.group(1))
        logs.add(re.search(r"^#SBATCH --output=(\S+)$", text, re.MULTILINE).group(1))
        ok("%s: uses the event manifest and matching case" % wrapper.name,
           "--manifest event_unbiased_manifest.json" in text
           and ("--case " + work.group(1)) in text)
        ok("%s: runs the analyzer after the MD" % wrapper.name,
           "analyze_event_unbiased.py" in text)
    ok("three distinct workdirs", workdirs == {"event_unbiased_s1", "event_unbiased_s2",
                                               "event_unbiased_s3"})
    ok("three distinct job names", len(job_names) == 3, str(sorted(job_names)))
    ok("three distinct log prefixes", len(logs) == 3)


def check_runner(cases):
    scratch = tempfile.mkdtemp(prefix="li2025_unbiased_selftest_")
    try:
        for case in cases:
            name = case["name"]
            workdir = Path(scratch) / name
            code, output = run_process([
                sys.executable, str(HERE / "run_cp2k.py"),
                "--manifest", "event_unbiased_manifest.json", "--case", name,
                "--workdir", str(workdir), "--workdir-root", scratch,
                "--ntasks", "64", "--dry-run"])
            payload = json.loads(output)
            ok("%s: dry-run stages two files and verifies hashes" % name,
               code == 0 and payload["dry_run"] and len(payload["staged"]) == 2)
            ok("%s: dry-run reports the seed and 3000-step budget" % name,
               payload["unbiased"]["seed"] == case["seed"]
               and payload["unbiased"]["steps"] == 3000
               and payload["unbiased"]["expected_frames"] == 301)
            ok("%s: staged input equals the generated input" % name,
               sha256(workdir / (name + ".inp"))
               == sha256(HERE / "inputs" / (name + ".inp")))
            ok("%s: no PLUMED file is staged" % name,
               not (workdir / (name + ".plumed")).exists())

        bulk = cases[0]
        workdir = Path(scratch) / "tamper"
        code, _ = run_process([
            sys.executable, str(HERE / "run_cp2k.py"),
            "--manifest", "event_unbiased_manifest.json", "--case", bulk["name"],
            "--workdir", str(workdir), "--workdir-root", scratch, "--dry-run"])
        ok("first dry-run into the tamper workdir succeeds", code == 0)
        xyz = workdir / "bulk_r1_event12270.xyz"
        xyz.write_text(xyz.read_text(encoding="utf-8") + "\n", encoding="utf-8")
        code, _ = run_process([
            sys.executable, str(HERE / "run_cp2k.py"),
            "--manifest", "event_unbiased_manifest.json", "--case", bulk["name"],
            "--workdir", str(workdir), "--workdir-root", scratch, "--dry-run"])
        ok("runner refuses a tampered event XYZ (exit 4)", code == 4)

        stale = Path(scratch) / "stale"
        stale.mkdir()
        (stale / (bulk["name"] + ".out")).write_text("old\n", encoding="utf-8")
        code, _ = run_process([
            sys.executable, str(HERE / "run_cp2k.py"),
            "--manifest", "event_unbiased_manifest.json", "--case", bulk["name"],
            "--workdir", str(stale), "--workdir-root", scratch, "--dry-run"])
        ok("runner refuses a stale output workdir (exit 7)", code == 7)

        code, _ = run_process([
            sys.executable, str(HERE / "run_cp2k.py"),
            "--manifest", "event_unbiased_manifest.json", "--case", bulk["name"],
            "--workdir", str(Path(scratch) / ".." / "escape"),
            "--workdir-root", scratch, "--dry-run"])
        ok("runner refuses '..' in the workdir (exit 5)", code == 5)

        code, _ = run_process([
            sys.executable, str(HERE / "run_cp2k.py"),
            "--manifest", "event_unbiased_manifest.json", "--case", bulk["name"],
            "--workdir", str(Path(scratch) / "nested" / "deep"),
            "--workdir-root", scratch, "--dry-run"])
        ok("runner refuses a nested workdir (exit 5)", code == 5)

        # End-to-end with fake CP2K output: complete (301 frames) and truncated.
        for label, last_step, frames, expected_code in (
                ("complete", 3000, 301, 0), ("short", 2999, 301, 6)):
            work = Path(scratch) / ("fake_" + label)

            def fake_cp2k(command, cwd=None, env=None, stdout=None, stderr=None,
                          last_step=last_step, frames=frames):
                Path(cwd, "event_unbiased_s1.out").write_text("\n".join([
                    " MD| Step number 1",
                    " SCF run converged in 8 steps",
                    " ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]:  -1735.10",
                    " MD| Step number %d" % last_step,
                    " SCF run converged in 9 steps",
                    " ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]:  -1735.11",
                    " PROGRAM ENDED AT   test",
                ]) + "\n", encoding="utf-8")
                body = []
                for index in range(frames):
                    step = index * 10
                    body.append("3\n i = %d, time = %d.000, E = -1.0\n"
                                "H 0 0 0\nH 0 0 1\nO 0 0 2\n" % (step, step))
                Path(cwd, "bulk_r1_unbiased-pos-1.xyz").write_text("".join(body),
                                                                   encoding="utf-8")
                Path(cwd, "bulk_r1_unbiased-1.ener").write_text(
                    "#     Step Nr.          Time[fs]        Kin.[a.u.]        "
                    "Temp[K]            Pot.[a.u.]        Cons Qty[a.u.]        UsedTime[s]\n"
                    "         0            0.000000         0.4       298.150000     "
                    "-1735.1       -1734.6       0.0\n"
                    "      %d          %d.000000         0.4       298.150000     "
                    "-1735.1       -1734.6       0.0\n" % (last_step, last_step),
                    encoding="utf-8")
                return 0

            original_call = runner.subprocess.call
            try:
                runner.subprocess.call = fake_cp2k
                code = runner.run("event_unbiased_s1", str(work), 64, False,
                                  "event_unbiased_manifest.json", None, scratch)
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
            else:
                ok("end-to-end truncated run exits 6 and records FAILED_STEP_RANGE",
                   code == 6 and record["status"] == "FAILED_STEP_RANGE"
                   and record["unbiased"]["step_check"] == "FAIL")
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def synthetic_frames(sc_values):
    base = []
    for index in range(96):
        x = (index % 12) * 1.1 + 0.5
        y = 10.0 + (index // 12) * 0.4
        z = 12.0
        base.append(("H", x, y, z))
        base.append(("H", x + 0.98, y, z))
        base.append(("O", x + 0.49, y + 0.85, z))
    cx, cy, cz = 2.0, 2.0, 2.0
    solute = [("C", cx, cy, cz), ("O", cx + 1.2, cy, cz),
              ("H", cx - 1.1, cy, cz), ("H", cx, cy + 1.1, cz)]
    frames = []
    for index, sc in enumerate(sc_values):
        sx = cx + sc
        atoms = list(base) + list(solute) + [
            ("S", sx, cy, cz),
            ("O", sx + 1.5, cy, cz),
            ("O", sx, cy + 1.5, cz),
            ("O", sx, cy, cz + 1.5),
        ]
        frames.append(atoms)
    return frames


def write_synthetic_trajectory(path, frames):
    lines = []
    for index, atoms in enumerate(frames):
        step = index * 10
        lines.append(str(len(atoms)))
        lines.append(" i = %d, time = %d.000, E = -1.0" % (step, step))
        for element, x, y, z in atoms:
            lines.append("%-2s %14.9f %14.9f %14.9f" % (element, x, y, z))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def check_analyzer():
    scratch = tempfile.mkdtemp(prefix="li2025_unbiased_analyzer_")
    try:
        values = [2.0] * 5 + [2.6] * 10 + [1.8] * 5
        trajectory = Path(scratch) / "synthetic.xyz"
        write_synthetic_trajectory(trajectory, synthetic_frames(values))
        out = Path(scratch) / "synthetic.json"
        code, output = run_process([
            sys.executable, str(HERE / "analyze_event_unbiased.py"),
            "--trajectory", str(trajectory), "--label", "synthetic",
            "--box", BULK_BOX, "--out", str(out)])
        ok("analyzer runs on a synthetic 20-frame trajectory", code == 0)
        data = json.loads(read_text(out))
        ok("analyzer counts 20 frames", data["frames"]["count"] == 20)
        summary = data["summary"]
        ok("S-C first/last/min/max are 2.0/1.8/1.8/2.6",
           abs(summary["sc"]["first"] - 2.0) < 1e-9
           and abs(summary["sc"]["last"] - 1.8) < 1e-9
           and abs(summary["sc"]["min"] - 1.8) < 1e-9
           and abs(summary["sc"]["max"] - 2.6) < 1e-9)
        ok("5 frames below 2 A with a longest run of 5",
           summary["sc"]["frames_below_threshold"] == 5
           and summary["sc"]["longest_below_threshold_frames"] == 5)
        ok("first sustained dissociation is frame 5 at 50 fs",
           summary["sc"]["first_sustained_dissociation"]["index"] == 5
           and summary["sc"]["first_sustained_dissociation"]["time_fs"] == 50.0)
        ok("C=O, C-H and S-O are reported from the solute geometry",
           abs(summary["co"]["first"] - 1.2) < 1e-6
           and abs(summary["ch1"]["first"] - 1.1) < 1e-6
           and abs(summary["ch2"]["first"] - 1.1) < 1e-6
           and abs(summary["so_mean"]["first"] - 1.5) < 1e-6)
        ok("nearest water H and its covalent O-H are reported",
           abs(data["per_frame"][0]["nearest_water_h_distance"] - 1.8) < 1e-6
           or data["per_frame"][0]["nearest_water_h_distance"] > 5.0)
        ok("the nearest-water O-H covalent distance is ~0.98 A",
           abs(data["per_frame"][0]["water_oh_distance"] - 0.981) < 1e-3,
           "%.4f" % data["per_frame"][0]["water_oh_distance"])
        ok("analyzer asserts no verdict",
           data["interpretation"]["verdict"] is None
           and "unbiased" in data["interpretation"]["what_this_is"].lower())
        ok("analyzer disclaims free energy and spontaneous formation",
           any("free energy" in item for item in data["interpretation"]["what_this_is_not"])
           and any("spontaneous" in item
                   for item in data["interpretation"]["what_this_is_not"]))
        tsv = out.with_name(out.stem + "_FRAMES.tsv")
        ok("analyzer writes a per-frame TSV",
           tsv.is_file() and len(read_text(tsv).splitlines()) == 21)
        ok("analyzer records the cell it was given",
           data["box_angstrom"] == [14.2296, 14.2296, 14.2296])
        code, _ = run_process([sys.executable, str(HERE / "analyze_event_unbiased.py"),
                               "--trajectory", str(trajectory),
                               "--out", str(Path(scratch) / "nobox.json")])
        ok("analyzer refuses to guess the cell when --box is omitted", code != 0)

        # No sustained dissociation when every frame stays above 2.5 A.
        flat = Path(scratch) / "flat.xyz"
        write_synthetic_trajectory(flat, synthetic_frames([1.8] * 20))
        flat_out = Path(scratch) / "flat.json"
        code, _ = run_process([sys.executable, str(HERE / "analyze_event_unbiased.py"),
                               "--trajectory", str(flat), "--box", BULK_BOX,
                               "--out", str(flat_out)])
        flat_data = json.loads(read_text(flat_out))
        ok("no sustained dissociation is reported when S-C never exceeds 2.5 A",
           code == 0 and flat_data["summary"]["sc"]["first_sustained_dissociation"] is None)

        # Temperature coverage.
        energy = Path(scratch) / "fake.ener"
        energy.write_text(
            "#     Step Nr.          Time[fs]        Kin.[a.u.]        Temp[K]            "
            "Pot.[a.u.]        Cons Qty[a.u.]        UsedTime[s]\n"
            "         0            0.000000         0.4       298.150000     -1.0       -0.6       0.0\n"
            "        10           10.000000         0.4       310.000000     -1.0       -0.6       0.0\n"
            "        20           20.000000         0.4       280.000000     -1.0       -0.6       0.0\n",
            encoding="utf-8")
        temp_out = Path(scratch) / "temp.json"
        code, _ = run_process([sys.executable, str(HERE / "analyze_event_unbiased.py"),
                               "--trajectory", str(trajectory), "--energy", str(energy),
                               "--box", BULK_BOX, "--out", str(temp_out)])
        temp = json.loads(read_text(temp_out))["temperature"]
        ok("temperature coverage is parsed from the CP2K .ener",
           code == 0 and temp["rows"] == 3 and abs(temp["mean_k"] - 296.05) < 1e-6
           and temp["min_k"] == 280.0 and temp["max_k"] == 310.0)

        # The seeded event XYZ itself must be analysable (no comma in its comment).
        single_out = Path(scratch) / "single.json"
        code, _ = run_process([sys.executable, str(HERE / "analyze_event_unbiased.py"),
                               "--trajectory", str(HERE / "inputs" / "bulk_r1_event12270.xyz"),
                               "--box", BULK_BOX, "--out", str(single_out)])
        single = json.loads(read_text(single_out)) if code == 0 else None
        ok("analyzer accepts the seeded event XYZ and reproduces S-C 1.778073 A",
           code == 0 and len(single["per_frame"]) == 1
           and abs(single["per_frame"][0]["sc"] - 1.778073) < 5e-4
           and abs(single["per_frame"][0]["nearest_water_h_distance"] - 1.478026) < 2e-3)

        # Bad atom order must be rejected.
        bad = Path(scratch) / "bad.xyz"
        lines = read_text(trajectory).splitlines()
        lines[2], lines[4] = lines[4], lines[2]
        bad.write_text("\n".join(lines) + "\n", encoding="utf-8")
        code, _ = run_process([sys.executable, str(HERE / "analyze_event_unbiased.py"),
                               "--trajectory", str(bad), "--box", BULK_BOX,
                               "--out", str(Path(scratch) / "bad.json")])
        ok("analyzer rejects a trajectory with the wrong atom order", code != 0)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def check_source_gates():
    scratch = tempfile.mkdtemp(prefix="li2025_unbiased_gates_")
    try:
        duplicate = Path(scratch) / "dup.xyz"
        frame = "1\n i = 100, time = 100.000, E = 0.0\nH 0 0 0\n"
        duplicate.write_text(frame + frame, encoding="utf-8")
        try:
            prep.find_event_frame(duplicate, 100.0)
            raise AssertionError("accepted two frames with the same time")
        except SystemExit:
            ok("duplicate frame times are rejected", True)

        missing = Path(scratch) / "missing.xyz"
        missing.write_text(frame.replace("100.000", "101.000"), encoding="utf-8")
        try:
            prep.find_event_frame(missing, 100.0)
            raise AssertionError("accepted a missing time")
        except SystemExit:
            ok("a missing frame time is rejected", True)

        small = Path(scratch) / "small.xyz"
        small.write_text("2\n i = 100, time = 100.000, E = 0.0\nH 0 0 0\nH 0 0 1\n",
                         encoding="utf-8")
        try:
            prep.find_event_frame(small, 100.0)
            raise AssertionError("accepted the wrong atom count")
        except SystemExit:
            ok("a frame with the wrong atom count is rejected", True)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def check_python36_and_shell():
    py_files = sorted(HERE.glob("*.py"))
    ok("the event scripts are present",
       {"prepare_event_unbiased.py", "analyze_event_unbiased.py",
        "selftest_event_unbiased.py"}.issubset(set(path.name for path in py_files)))
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


def check_scope():
    readme = read_text(HERE / "README.md").lower()
    for token in ("unbiased", "geoopt", "slot", "verbatim", "not a free energy"):
        ok("README documents %r" % token, token in readme)
    leftovers = [path.name for path in HERE.rglob("*")
                 if path.suffix.lower() in (".hills", ".colvar")]
    ok("no HILLS/COLVAR files in the campaign directory", not leftovers, str(leftovers))


def main():
    manifest = json.loads(read_text(MANIFEST))
    geometry = check_source(manifest)
    cases = check_manifest(manifest)
    check_inputs(cases, geometry)
    check_wrappers()
    check_runner(cases)
    check_analyzer()
    check_source_gates()
    check_python36_and_shell()
    check_scope()
    print("SELFTEST_EVENT_UNBIASED_OK: %d checks passed (source, inputs, wrappers, "
          "runner, analyzer, gates, py36, scope)" % CHECKS["count"])
    return 0


if __name__ == "__main__":
    sys.exit(main())

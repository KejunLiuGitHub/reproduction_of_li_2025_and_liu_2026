#!/usr/bin/env python3
"""Self-test for the Li 2025 direct second segment (seg2) continuation.

Uses the real pulled-back parents for integration checks and synthetic strings
for the hard gates, without touching the cluster and without running CP2K.

Covered:
  * parent state: COMPLETED run_record, job ids, &MD STEP_START_VAL 10000 (with
    the &METADYN STEP_START_VAL 0 decoy scoped out), 296 COORD/VELOCITY rows,
    HILLS/COLVAR baselines, energy/trajectory ending at 10000,
  * generated CP2K input: one &EXT_RESTART block immediately after &END GLOBAL,
    RESTART_DEFAULT .FALSE., POS/VEL/COUNTERS .TRUE., bare restart basename,
    STEPS unchanged and every scientific line byte-identical,
  * generated PLUMED input: one RESTART action as the first action after UNITS,
    same CV/SIGMA/PACE/HEIGHT/wall and same HILLS/COLVAR file names,
  * staging plan: nine files, non-colliding targets, all SHA256 re-derived,
  * runner: dry-run staging for both cases, runtime STEP_START_VAL gate,
    hash-mismatch refusal (exit 4), stale-output refusal (exit 7) and unsafe
    workdir refusal (exit 5, including '..', nesting and workdir == parent),
  * Python 3.6 syntax/API and `bash -n` for every script in the directory.

Run: python3 selftest_segment2.py
Exit 0 and a final SELFTEST_OK line on success, 1 on any failure.
"""
import sys

sys.dont_write_bytecode = True

import ast
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

import prepare_segment2 as prep
import run_cp2k as runner

PARENT_ROOT = ROOT / "results" / "li2025_direct"
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


def run_runner(args):
    process = subprocess.Popen(
        [sys.executable, str(HERE / "run_cp2k.py")] + args,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = process.communicate()[0].decode("utf-8", "replace")
    return process.returncode, output


def load_manifest():
    return json.loads(read_text(HERE / "seg2_manifest.json"))


def check_manifest(manifest):
    ok("seg2 manifest schema", manifest.get("schema") == "li2025_direct_seg2/v1")
    cases = manifest.get("cases", [])
    ok("seg2 manifest has exactly two cases",
       [case["name"] for case in cases] == ["bulk_r1_seg2", "interface_r1_seg2"])

    expected_jobs = {"bulk_r1_seg2": "8562691", "interface_r1_seg2": "8562695"}
    for case in cases:
        name = case["name"]
        segment = case["segment"]
        ok("%s: parent job id %s" % (name, expected_jobs[name]),
           case["parent_job_id"] == expected_jobs[name])
        ok("%s: segment is steps 10001-20000 (cumulative 1-20000)" % name,
           segment["start_step"] == 10000
           and segment["steps"] == 10000
           and segment["expected_first_step"] == 10001
           and segment["expected_end_step"] == 20000
           and segment["segment_step_range"] == [10001, 20000]
           and segment["cumulative_step_range"] == [1, 20000]
           and segment["parent_md_step_last"] == 10000)
        ok("%s: still no FES" % name, case["fes"]["status"] == "not_computed")
        ok("%s: PLUMED RESTART was added because the parent had none" % name,
           case["plumed_restart_added"] is True)

        evidence = case["parent_evidence"]
        ok("%s: parent run_record is COMPLETED with 10000 steps" % name,
           evidence["run_record"]["status"] == "COMPLETED"
           and evidence["run_record"]["program_ended"] is True
           and evidence["run_record"]["md_step_last"] == 10000
           and evidence["run_record"]["scf_not_converged_runs"] == 0)
        restart = evidence["restart"]
        ok("%s: restart &MD STEP_START_VAL is 10000 at ~10000 fs" % name,
           restart["md_step_start_val"] == 10000
           and abs(restart["md_time_start_val_fs"] - 10000.0) < 1.0e-3
           and restart["md_steps_field"] == 10000)
        ok("%s: restart carries 296 COORD and 296 VELOCITY rows" % name,
           restart["coord_atom_lines"] == 296 and restart["velocity_atom_lines"] == 296)
        ok("%s: restart also contains the &METADYN STEP_START_VAL 0 decoy" % name,
           set(restart["all_step_start_val_values"]) == {"0", "10000"},
           str(restart["all_step_start_val_values"]))
        ok("%s: parent energy and trajectory end at step 10000" % name,
           evidence["energy"]["last_step"] == 10000
           and evidence["trajectory"]["last_step"] == 10000)
        baseline = case["metad_baseline"]
        ok("%s: parent HILLS/COLVAR baselines are non-empty and below 10000 fs" % name,
           baseline["hills_rows"] > 0 and baseline["colvar_rows"] > 0
           and baseline["hills_last_time_fs"] < 10000
           and baseline["colvar_last_time_fs"] < 10000)
        ok("%s: parent workdir is the fix1 directory" % name,
           case["parent_workdir_name"] == case["parent_case"] + "_fix1")
    return cases


def check_generated_inputs(case):
    name = case["name"]
    parent_name = case["parent_case"]
    parent_dir = PARENT_ROOT / case["parent_workdir_name"]
    parent_inp = read_text(parent_dir / (parent_name + ".inp"))
    parent_plumed = read_text(parent_dir / (parent_name + ".plumed"))
    new_inp = read_text(HERE / "inputs" / (name + ".inp"))
    new_plumed = read_text(HERE / "inputs" / (name + ".plumed"))

    restart_target = next(entry["target"] for entry in case["stage"]
                          if entry["role"] == "parent_restart")
    block = prep.EXT_RESTART_TEMPLATE.format(restart=restart_target)
    ok("%s: input is the parent plus exactly one &EXT_RESTART block" % name,
       new_inp.replace(block, "", 1) == parent_inp and new_inp.count("&EXT_RESTART") == 1)
    lines = new_inp.splitlines()
    index = lines.index("&END GLOBAL")
    ok("%s: &EXT_RESTART sits immediately after &END GLOBAL" % name,
       lines[index + 1] == "&EXT_RESTART"
       and lines[index + 2] == "  RESTART_FILE_NAME " + restart_target
       and lines.index("&END EXT_RESTART") < lines.index("&FORCE_EVAL"))
    ok("%s: EXT_RESTART fields are the required ones" % name,
       "RESTART_DEFAULT .FALSE." in block
       and "RESTART_POS .TRUE." in block
       and "RESTART_VEL .TRUE." in block
       and "RESTART_COUNTERS .TRUE." in block)
    ok("%s: restart reference is a bare basename inside the workdir" % name,
       "/" not in restart_target and "\\" not in restart_target
       and ".." not in restart_target
       and ("RESTART_FILE_NAME " + restart_target) in block)
    ok("%s: STEPS is still 10000 (per-segment budget)" % name,
       prep.steps_line_value(new_inp) == 10000 and new_inp.count("STEPS 10000") == 1)
    for token in ("CHARGE -2", "MULTIPLICITY 1", "&XC_FUNCTIONAL BLYP",
                  "TYPE DFTD3(BJ)", "REFERENCE_FUNCTIONAL BLYP", "TIMESTEP 1.0",
                  "TEMPERATURE 298.15", "TIMECON 100", "MD 10", "MD 100"):
        ok("%s: scientific line %r unchanged" % (name, token), token in new_inp)

    units_line = next(line for line in parent_plumed.splitlines(True)
                      if line.upper().startswith("UNITS "))
    ok("%s: PLUMED is the parent plus exactly one RESTART action" % name,
       new_plumed == parent_plumed.replace(units_line, units_line + "RESTART\n", 1)
       and sum(1 for line in new_plumed.splitlines()
               if line.strip().upper() == "RESTART") == 1)
    plumed_lines = new_plumed.splitlines()
    units_index = [i for i, line in enumerate(plumed_lines)
                   if line.upper().startswith("UNITS ")][0]
    ok("%s: RESTART is the first action after UNITS" % name,
       plumed_lines[units_index + 1].strip() == "RESTART")
    for token in ("DISTANCE ATOMS=293,289", "SIGMA=0.1", "PACE=30",
                  "FILE=HILLS_" + parent_name, "FILE=COLVAR_" + parent_name,
                  "UPPER_WALLS ARG=d_sc AT=4.0 KAPPA=200.0"):
        ok("%s: PLUMED term %r unchanged" % (name, token), token in new_plumed)
    ok("%s: HILLS/COLVAR names are unchanged so PLUMED appends" % name,
       case["metad_files"]["hills"] == "HILLS_" + parent_name
       and case["metad_files"]["colvar"] == "COLVAR_" + parent_name)


def check_stage_plan(case):
    name = case["name"]
    entries = case["stage"]
    ok("%s: nine staged files" % name, len(entries) == 9)
    roles = [entry["role"] for entry in entries]
    ok("%s: required roles are present" % name,
       roles == ["new_input", "new_plumed", "parent_restart", "parent_xyz",
                 "parent_input", "parent_plumed", "parent_run_record",
                 "parent_hills", "parent_colvar"])
    parent_dir = PARENT_ROOT / case["parent_workdir_name"]
    targets = [entry["target"] for entry in entries]
    ok("%s: staged targets are unique" % name, len(set(targets)) == len(targets))
    for entry in entries:
        base = parent_dir if entry["role"].startswith("parent") else HERE / "inputs"
        source = base / entry["source"]
        ok("%s: %s sha256 matches" % (name, entry["target"]),
           source.is_file() and sha256(source) == entry["sha256"])
    outputs = case["outputs"]
    ok("%s: no staged file collides with a CP2K output" % name,
       outputs["restart"] not in targets and outputs["trajectory"] not in targets
       and outputs["energy"] not in targets and outputs["output"] not in targets)
    ok("%s: copied parent restart has a non-colliding name" % name,
       next(entry["target"] for entry in entries
            if entry["role"] == "parent_restart") == name + "_parent.restart")
    ok("%s: CP2K project and output names are the parent's" % name,
       case["cp2k_project"] == case["parent_case"]
       and outputs["trajectory"] == case["parent_case"] + "-pos-1.xyz"
       and outputs["energy"] == case["parent_case"] + "-1.ener"
       and outputs["restart"] == case["parent_case"] + "-1.restart")


def check_runner(cases):
    scratch = tempfile.mkdtemp(prefix="li2025_seg2_selftest_")
    try:
        for case in cases:
            name = case["name"]
            parent_dir = PARENT_ROOT / case["parent_workdir_name"]
            workdir = Path(scratch) / name
            code, output = run_runner([
                "--manifest", "seg2_manifest.json", "--case", name,
                "--parent-dir", str(parent_dir), "--workdir", str(workdir),
                "--workdir-root", scratch, "--ntasks", "64", "--dry-run"])
            payload = json.loads(output)
            ok("%s: runner dry-run stages and hash-checks nine files" % name,
               code == 0 and payload["dry_run"] and len(payload["staged"]) == 9)
            ok("%s: command launches 64-rank cp2k.psmp on the seg2 input" % name,
               payload["command"][1:3] == ["-np", "64"]
               and payload["command"][3] == runner.CP2K
               and payload["command"][4:6] == ["-i", name + ".inp"]
               and payload["command"][6:8] == ["-o", name + ".out"])
            ok("%s: staged restart still starts at 10000" % name,
               runner.restart_md_start_step(workdir / (name + "_parent.restart")) == 10000)
            staged_plumed = read_text(workdir / (case["parent_case"] + ".plumed"))
            ok("%s: staged PLUMED carries RESTART" % name,
               any(line.strip().upper() == "RESTART" for line in staged_plumed.splitlines()))
            ok("%s: staged HILLS/COLVAR are the parent files" % name,
               sha256(workdir / case["metad_files"]["hills"])
               == sha256(parent_dir / case["metad_files"]["hills"])
               and sha256(workdir / case["metad_files"]["colvar"])
               == sha256(parent_dir / case["metad_files"]["colvar"]))

        bulk = cases[0]
        parent_dir = PARENT_ROOT / bulk["parent_workdir_name"]

        # Hash mismatch: tamper with the staged HILLS file.
        workdir = Path(scratch) / "tamper"
        code, _ = run_runner([
            "--manifest", "seg2_manifest.json", "--case", bulk["name"],
            "--parent-dir", str(parent_dir), "--workdir", str(workdir),
            "--workdir-root", scratch, "--dry-run"])
        ok("dry-run into a fresh workdir succeeds", code == 0)
        hills = workdir / bulk["metad_files"]["hills"]
        hills.write_text(hills.read_text(encoding="utf-8") + "\n", encoding="utf-8")
        code, _ = run_runner([
            "--manifest", "seg2_manifest.json", "--case", bulk["name"],
            "--parent-dir", str(parent_dir), "--workdir", str(workdir),
            "--workdir-root", scratch, "--dry-run"])
        ok("runner refuses a tampered HILLS file (exit 4)", code == 4)

        # Stale output: a previous CP2K output must not be overwritten.
        stale = Path(scratch) / "stale"
        stale.mkdir()
        (stale / (bulk["name"] + ".out")).write_text("old run\n", encoding="utf-8")
        code, _ = run_runner([
            "--manifest", "seg2_manifest.json", "--case", bulk["name"],
            "--parent-dir", str(parent_dir), "--workdir", str(stale),
            "--workdir-root", scratch, "--dry-run"])
        ok("runner refuses a workdir with stale output (exit 7)", code == 7)

        # Unsafe paths.
        code, _ = run_runner([
            "--manifest", "seg2_manifest.json", "--case", bulk["name"],
            "--parent-dir", str(parent_dir), "--workdir",
            str(Path(scratch) / ".." / "escape"), "--workdir-root", scratch, "--dry-run"])
        ok("runner refuses '..' in the workdir (exit 5)", code == 5)

        code, _ = run_runner([
            "--manifest", "seg2_manifest.json", "--case", bulk["name"],
            "--parent-dir", str(parent_dir), "--workdir",
            str(Path(scratch) / "nested" / "deep"), "--workdir-root", scratch, "--dry-run"])
        ok("runner refuses a nested workdir (exit 5)", code == 5)

        code, _ = run_runner([
            "--manifest", "seg2_manifest.json", "--case", bulk["name"],
            "--parent-dir", str(parent_dir), "--workdir", str(parent_dir),
            "--workdir-root", str(PARENT_ROOT), "--dry-run"])
        ok("runner refuses to use the parent workdir as the segment workdir (exit 5)",
           code == 5)

        # Runtime STEP_START_VAL gate: a manifest claiming a different start must fail.
        tampered_manifest = json.loads(read_text(HERE / "seg2_manifest.json"))
        for entry in tampered_manifest["cases"]:
            if entry["name"] == bulk["name"]:
                entry["segment"]["start_step"] = 9999
        # The runner only accepts registered manifest file names, so the tampered
        # copy keeps the registered name inside its own directory.
        tampered_dir = Path(scratch) / "tampered_manifest"
        tampered_dir.mkdir()
        manifest_path = tampered_dir / "seg2_manifest.json"
        manifest_path.write_text(json.dumps(tampered_manifest, indent=2), encoding="utf-8")
        code, _ = run_runner([
            "--manifest", str(manifest_path), "--case", bulk["name"],
            "--parent-dir", str(parent_dir), "--workdir", str(Path(scratch) / "badstart"),
            "--workdir-root", scratch, "--dry-run"])
        ok("runner refuses a manifest whose start step disagrees with the checkpoint (exit 4)",
           code == 4)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def check_segment_audit(manifest):
    case = manifest["cases"][0]
    parsed = {"md_step_first": 10001, "md_step_last": 20000}
    audit = runner.audit_segment(case, parsed, 1000, [10010, 20000], 20000,
                                 (666, 19979.999), (2000, 19999.999), 10000)
    ok("audit: a complete segment passes step, energy, trajectory and append checks",
       audit["step_check"] == "PASS" and audit["energy_check"] == "PASS"
       and audit["trajectory"]["check"] == "PASS"
       and audit["hills"]["check"] == "PASS" and audit["colvar"]["check"] == "PASS"
       and audit["hills"]["rows_segment"] == 333
       and audit["colvar"]["rows_segment"] == 1000
       and audit["trajectory"]["no_overlap_with_parent"] is True)
    bad = runner.audit_segment(case, {"md_step_first": 10001, "md_step_last": 19999},
                               1, [19999], 19999, (333, 9989.999), (1000, 9999.999), 10000)
    ok("audit: a truncated segment fails the step, trajectory and append checks",
       bad["step_check"] == "FAIL" and bad["energy_check"] == "FAIL"
       and bad["trajectory"]["check"] == "FAIL"
       and bad["hills"]["check"] == "FAIL" and bad["colvar"]["check"] == "FAIL")

    scratch = tempfile.mkdtemp(prefix="li2025_seg2_e2e_")
    workdir = Path(scratch) / case["name"]
    parent_dir = PARENT_ROOT / case["parent_workdir_name"]
    original_call = runner.subprocess.call

    def fake_cp2k(command, cwd=None, env=None, stdout=None, stderr=None):
        out = Path(cwd) / "bulk_r1_seg2.out"
        out.write_text("\n".join([
            " MD| Step number 10001",
            " SCF run converged in 8 steps",
            " ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]:  -1735.100000",
            " MD| Step number 20000",
            " SCF run converged in 9 steps",
            " ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]:  -1735.110000",
            " **** ******  **  PROGRAM ENDED AT   test",
        ]) + "\n", encoding="utf-8")
        trajectory = Path(cwd) / "bulk_r1-pos-1.xyz"
        frame = "3\n i = {step}, time = {time}, E = -1.0\nH 0 0 0\nH 0 0 1\nO 0 0 2\n"
        trajectory.write_text(frame.format(step=10010, time="10010.000")
                              + frame.format(step=20000, time="20000.000"),
                              encoding="utf-8")
        (Path(cwd) / "bulk_r1-1.ener").write_text(
            "# header\n19999 19999.0 0.4 300.0 -1735.1 -1734.6 6.8\n"
            "20000 20000.0 0.4 300.0 -1735.11 -1734.6 6.8\n", encoding="utf-8")
        with open(str(Path(cwd) / "HILLS_bulk_r1"), "a") as handle:
            handle.write("19979.999 3.20 0.1 0.0627509 -1\n")
        with open(str(Path(cwd) / "COLVAR_bulk_r1"), "a") as handle:
            handle.write("19999.999 3.20 1.000000\n")
        return 0

    try:
        runner.subprocess.call = fake_cp2k
        code = runner.run(case["name"], str(workdir), 64, False, "seg2_manifest.json",
                          str(parent_dir), scratch)
    finally:
        runner.subprocess.call = original_call
    try:
        record = json.loads(read_text(workdir / "run_record.json"))
        ok("end-to-end: a complete fake segment exits 0 and records COMPLETED",
           code == 0 and record["status"] == "COMPLETED"
           and record["segment"]["step_check"] == "PASS"
           and record["segment"]["energy_check"] == "PASS"
           and record["segment"]["trajectory"]["check"] == "PASS"
           and record["segment"]["hills"]["check"] == "PASS"
           and record["segment"]["colvar"]["check"] == "PASS")
        ok("end-to-end: run_record carries parent job, ranges and no FES",
           record["segment"]["parent_job_id"] == case["parent_job_id"]
           and record["segment"]["start_step"] == 10000
           and record["segment"]["cumulative_step_range"] == [1, 20000]
           and record["fes"]["status"] == "not_computed"
           and record["scf"]["not_converged_runs"] == 0)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def check_synthetic_gates():
    synthetic_restart = "\n".join([
        " # restart",
        " &MOTION",
        "   &MD",
        "     STEPS  10000",
        "     STEP_START_VAL  10000",
        "     TIME_START_VAL  9.9999999999989159E+03",
        "   &END MD",
        "   &FREE_ENERGY",
        "     &METADYN",
        "       STEP_START_VAL  0",
        "     &END METADYN",
        "   &END FREE_ENERGY",
        " &END MOTION",
        "",
    ])
    counters = prep.md_block_counters(synthetic_restart)
    ok("synthetic restart: &MD counters ignore the &METADYN STEP_START_VAL 0 decoy",
       counters["STEP_START_VAL"] == "10000"
       and abs(float(counters["TIME_START_VAL"]) - 10000.0) < 1.0e-3)
    ok("synthetic restart: scoped atom counter works for COORD and VELOCITY",
       prep.section_atom_lines("&SUBSYS\n&COORD\nH 0 0 0\nO 1 0 0\n&END COORD\n"
                               "&VELOCITY\n1.0 2.0 3.0\n4.0 5.0 6.0\n&END VELOCITY\n&END SUBSYS\n",
                               "COORD") == 2
       and prep.section_atom_lines("&SUBSYS\n&COORD\nH 0 0 0\nO 1 0 0\n&END COORD\n"
                                   "&VELOCITY\n1.0 2.0 3.0\n4.0 5.0 6.0\n&END VELOCITY\n"
                                   "&END SUBSYS\n", "VELOCITY") == 2)

    inp = "&GLOBAL\n  PROJECT x\n&END GLOBAL\n&FORCE_EVAL\n&END FORCE_EVAL\n"
    inserted = prep.insert_ext_restart(inp, "x_parent.restart")
    lines = inserted.splitlines()
    ok("synthetic input: &EXT_RESTART is inserted right after &END GLOBAL",
       lines[2] == "&END GLOBAL" and lines[3] == "&EXT_RESTART"
       and inserted.index("&END EXT_RESTART") < inserted.index("&FORCE_EVAL"))
    for bad in ("../x.restart", "dir/x.restart", "x/../y.restart"):
        try:
            prep.insert_ext_restart(inp, bad)
            raise AssertionError("accepted unsafe restart name %r" % bad)
        except ValueError:
            pass
    ok("synthetic input: restart names with directories or '..' are rejected", True)
    try:
        prep.insert_ext_restart(inserted, "y.restart")
        raise AssertionError("accepted a second &EXT_RESTART")
    except ValueError:
        ok("synthetic input: a second &EXT_RESTART is rejected", True)
    try:
        prep.insert_ext_restart("&GLOBAL\n&END FORCE_EVAL\n", "x.restart")
        raise AssertionError("accepted an input without &END GLOBAL")
    except ValueError:
        ok("synthetic input: missing &END GLOBAL is rejected", True)

    plumed = "UNITS LENGTH=A ENERGY=kcal/mol TIME=fs\nd: DISTANCE ATOMS=1,2\n"
    with_restart, added = prep.insert_plumed_restart(plumed)
    ok("synthetic PLUMED: RESTART is inserted after UNITS", added
       and with_restart.splitlines()[1] == "RESTART")
    unchanged, added_again = prep.insert_plumed_restart(with_restart)
    ok("synthetic PLUMED: an existing RESTART is kept, not duplicated",
       not added_again and unchanged == with_restart)
    no_units, _ = prep.insert_plumed_restart("d: DISTANCE ATOMS=1,2\n")
    ok("synthetic PLUMED: without UNITS, RESTART goes on line 1",
       no_units.splitlines()[0] == "RESTART")

    try:
        prep.steps_line_value("&MD\n  STEPS 10000\n  STEPS 10000\n&END MD\n")
        raise AssertionError("accepted duplicate STEPS lines")
    except SystemExit:
        ok("synthetic input: duplicate STEPS lines are rejected", True)


def check_wrappers():
    wrappers = {
        "bulk_r1_seg2": HERE / "submit_bulk_r1_seg2.slurm",
        "interface_r1_seg2": HERE / "submit_interface_r1_seg2.slurm",
    }
    workdirs = set()
    for name, path in wrappers.items():
        text = read_text(path)
        ok("%s: wrapper exists and passes bash -n" % name,
           path.is_file() and subprocess.call(["bash", "-n", str(path)]) == 0)
        ok("%s: 64 MPI ranks, 32G, 48 h, mail END,FAIL" % name,
           "#SBATCH --ntasks=64" in text and "#SBATCH --mem=32G" in text
           and "#SBATCH --time=48:00:00" in text
           and "#SBATCH --mail-type=END,FAIL" in text)
        ok("%s: no GPU request and one thread per rank" % name,
           "--gres" not in text and "--cpus-per-task" not in text
           and "OMP_NUM_THREADS=1" in text)
        ok("%s: wrapper uses the seg2 manifest and case" % name,
           "--manifest seg2_manifest.json" in text and "--case %s" % name in text)
        match = re.search(r"^WORK=\$BASE/(\S+)$", text, re.MULTILINE)
        ok("%s: independent workdir" % name, match is not None and match.group(1) == name)
        if match:
            workdirs.add(match.group(1))
        ok("%s: parent COMPLETED guard is present" % name,
           "run_record.json" in text and "COMPLETED" in text)
    ok("the two seg2 wrappers use distinct workdirs",
       workdirs == {"bulk_r1_seg2", "interface_r1_seg2"})


def check_python36_and_shell():
    py_files = sorted(HERE.glob("*.py"))
    names = set(path.name for path in py_files)
    ok("all expected python entry points exist",
       {"prepare_li2025_direct.py", "prepare_segment2.py", "run_cp2k.py",
        "selftest_li2025_direct.py", "selftest_segment2.py"}.issubset(names),
       str(sorted(names)))
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
        subprocess_pattern = (r"subprocess\.(run|Popen|call|check_call|check_output)"
                              r"\([^)]*?\b(" + text_keyword + "|" + capture_keyword + r")\s*=")
        ok("%s has no subprocess text/capture keywords" % path.name,
           re.search(subprocess_pattern, source, re.DOTALL) is None)
        ok("%s has no walrus operator" % path.name, (":" + "=") not in source)
    for path in sorted(HERE.glob("*.slurm")):
        ok("%s passes bash -n" % path.name, subprocess.call(["bash", "-n", str(path)]) == 0)


def check_scope(manifest):
    readme = read_text(HERE / "README.md").lower()
    for token in ("seg2", "ext_restart", "restart_default .false.", "stale",
                  "not a converged fes", "step_start_val"):
        ok("README documents %r" % token, token in readme)
    policy = manifest.get("continuation_policy", {})
    ok("manifest states STEPS is a per-segment budget",
       "per-segment budget" in policy.get("steps_semantics", ""))
    ok("manifest states PLUMED appends instead of restarting from zero",
       "appends" in policy.get("plumed", ""))
    leftovers = [path.name for path in HERE.rglob("*")
                 if path.suffix.lower() in (".hills", ".colvar")
                 or path.name.startswith(("HILLS_", "COLVAR_"))]
    ok("no HILLS/COLVAR artefacts in the campaign directory", not leftovers, str(leftovers))


def main():
    manifest = load_manifest()
    cases = check_manifest(manifest)
    for case in cases:
        check_generated_inputs(case)
        check_stage_plan(case)
    check_runner(cases)
    check_segment_audit(manifest)
    check_synthetic_gates()
    check_wrappers()
    check_python36_and_shell()
    check_scope(manifest)
    print("SELFTEST_SEG2_OK: %d checks passed (parent state, inputs, staging, runner, "
          "synthetic gates, wrappers, py36, scope)" % CHECKS["count"])
    return 0


if __name__ == "__main__":
    sys.exit(main())

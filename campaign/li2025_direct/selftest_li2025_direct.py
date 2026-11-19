#!/usr/bin/env python3
"""Self-test for campaign/li2025_direct (Li 2025 Figure 6b pathway 3, first replicas).

Verifies, without touching the cluster and without running CP2K:

  * the workbook MD5 and that both exported blocks still match the sheet,
  * composition, atom order, S/C 1-based indices, S-C distances and contacts,
  * verbatim coordinate fidelity of the two XYZ files,
  * CP2K settings (box, BLYP + DFTD3(BJ)/BLYP, basis/potentials, cutoff,
    charge -2 / multiplicity 1, NVT/CSVR, 1 fs, 10 ps, SCF settings),
  * PLUMED settings (single d(S-C), units, SIGMA, PACE, HEIGHT conversion,
    absence of BIASFACTOR, UPPER_WALLS values and the R-snapshot wall margin),
  * SLURM wrappers (64 MPI ranks, 32G, 48 h, mail END/FAIL, no --gres),
  * runner dry-run staging, hash verification and tamper rejection,
  * Python 3.6 syntax/API compatibility and `bash -n` on the wrappers,
  * the scope disclaimers (no FES, not Li's three replicas, not Liu's five-CV).

Run: python3 selftest_li2025_direct.py
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
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import prepare_li2025_direct as prep
import run_cp2k as runner

CHECKS = {"count": 0}


def ok(label, condition, detail=""):
    CHECKS["count"] += 1
    if not condition:
        raise AssertionError("%s%s" % (label, (": " + detail) if detail else ""))
    print("ok %2d - %s" % (CHECKS["count"], label))


def read_xyz(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    count = int(lines[0])
    if len(lines) != count + 2:
        raise AssertionError("%s: %d lines for %d atoms" % (path, len(lines), count))
    elements = []
    coords = []
    for line in lines[2:]:
        parts = line.split()
        elements.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return elements, coords


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class StubWorksheet(object):
    def __init__(self, rows):
        self._rows = rows

    def iter_rows(self, min_row, max_row, values_only):
        return iter(self._rows)


def check_source_and_export():
    ok("workbook MD5 is the pinned value",
       prep.md5(prep.SOURCE) == prep.SOURCE_MD5, prep.md5(prep.SOURCE))

    from openpyxl import load_workbook
    workbook = load_workbook(str(prep.SOURCE), read_only=True, data_only=True)
    sheet = workbook[prep.SHEET]
    header = [cell.value for cell in next(sheet.iter_rows(min_row=1, max_row=1))]
    ok("sheet Figure 6b has the six expected block headers",
       [header[index] for index in (0, 5, 10, 15, 20, 25)]
       == ["R_bulk", "TS_bulk", "P_bulk", "R_inter", "TS_inter", "P_inter"])

    exported = {}
    for label, column in prep.BLOCKS:
        elements, coords = prep.read_block(sheet, column)
        box = {"R_bulk": (14.2296,) * 3, "R_inter": (15.0, 15.0, 30.0)}[label]
        metrics = prep.verify_structure(elements, coords, box, label)
        exported[label] = (elements, coords, metrics)
        ok("%s: 296 atoms H194 O100 C1 S1 in H,H,O x96 + HCHO + SO3 order" % label,
           len(elements) == 296
           and dict(Counter(elements)) == {"H": 194, "O": 100, "C": 1, "S": 1}
           and elements == prep.EXPECTED_ORDER)
        ok("%s: S is atom 293 and C is atom 289 (1-based)" % label,
           elements.index("S") == 292 and elements.index("C") == 288)
        ok("%s: no exact duplicate atoms (closest pair %.4f A)"
           % (label, metrics["min_atom_pair_distance_angstrom"]),
           metrics["min_atom_pair_distance_angstrom"] > 1.0e-6)
    workbook.close()

    ok("R_bulk S-C = 3.3191 A", abs(exported["R_bulk"][2]["sc_distance_angstrom"] - 3.3191) < 5e-4,
       "%.4f" % exported["R_bulk"][2]["sc_distance_angstrom"])
    ok("R_inter S-C = 3.9317 A",
       abs(exported["R_inter"][2]["sc_distance_angstrom"] - 3.9317) < 5e-4,
       "%.4f" % exported["R_inter"][2]["sc_distance_angstrom"])
    ok("R_bulk min intermolecular contact = 1.4623 A",
       abs(exported["R_bulk"][2]["min_intermolecular_contact_angstrom"] - 1.4623) < 1e-3,
       "%.4f" % exported["R_bulk"][2]["min_intermolecular_contact_angstrom"])
    ok("R_inter min intermolecular contact = 1.5490 A",
       abs(exported["R_inter"][2]["min_intermolecular_contact_angstrom"] - 1.5490) < 1e-3,
       "%.4f" % exported["R_inter"][2]["min_intermolecular_contact_angstrom"])

    # Hard-gate negative controls: duplicates and unparseable rows must be rejected.
    elements, coords, _ = exported["R_bulk"]
    duplicated = [list(point) for point in coords]
    duplicated[1] = list(duplicated[0])
    try:
        prep.verify_structure(elements, duplicated, (14.2296,) * 3, "dup")
        raise AssertionError("verify_structure accepted an exact duplicate atom")
    except ValueError:
        ok("exact duplicate atoms are rejected", True)

    rows = [("H", 0.0, 0.0, 0.0)] * 296
    broken = list(rows)
    broken[5] = (None, 0.0, 0.0, 0.0)
    try:
        prep.read_block(StubWorksheet(broken), 1)
        raise AssertionError("read_block accepted a missing element")
    except ValueError:
        ok("unparseable rows are rejected", True)
    return exported


def check_manifest(exported):
    manifest = json.loads((HERE / "manifest.json").read_text(encoding="utf-8"))
    ok("manifest pins the workbook MD5", manifest["source"]["workbook_md5"] == prep.SOURCE_MD5)
    ok("manifest records sheet and row range",
       manifest["source"]["sheet"] == "Figure 6b"
       and manifest["source"]["rows"] == "2:297"
       and manifest["source"]["atom_count_per_block"] == 296)
    ok("manifest records the exported cell ranges",
       manifest["source"]["blocks_exported"]["R_bulk"]["cell"] == "A2:D297"
       and manifest["source"]["blocks_exported"]["R_inter"]["cell"] == "P2:S297")
    ok("manifest system is q=-2, multiplicity 1, no counterions",
       manifest["system"]["net_charge"] == -2
       and manifest["system"]["multiplicity"] == 1
       and manifest["system"]["counterions"].startswith("none"))
    ok("manifest says the q=-2 model must not be mixed with the +2K model",
       "not the project's +2 K+ neutral model" in manifest["system"]["charge_model_note"])

    cases = {record["name"]: record for record in manifest["cases"]}
    ok("manifest has exactly the two R cases", sorted(cases) == ["bulk_r1", "interface_r1"])
    ok("no FES is claimed", manifest["fes"]["status"] == "not_computed")
    scope = " ".join(manifest["scope"]["what_this_is_not"]).lower()
    ok("scope disclaims Li's three replicas and Liu's five-distance CV",
       "three-replica" in scope and "five-distance" in scope)

    for name in ("bulk_r1", "interface_r1"):
        record = cases[name]
        point = record["point"]
        source_elements, source_coords, metrics = exported[point]
        elements, coords = read_xyz(HERE / record["xyz"])
        ok("%s: XYZ has 296 atoms in the source order" % name,
           elements == source_elements)
        max_delta = max(abs(coords[i][k] - source_coords[i][k])
                        for i in range(296) for k in range(3))
        ok("%s: coordinates are verbatim (max delta %.2e A)" % (name, max_delta),
           max_delta <= 1.0e-6)
        ok("%s: manifest SHA256 matches the three input files" % name,
           sha256(HERE / record["xyz"]) == record["xyz_sha256"]
           and sha256(HERE / record["inp"]) == record["inp_sha256"]
           and sha256(HERE / record["plumed"]) == record["plumed_sha256"])
        ok("%s: manifest S-C matches the re-extracted value" % name,
           abs(record["sc_distance_angstrom"] - metrics["sc_distance_angstrom"]) < 1.0e-9)
        ok("%s: S-C is inside the 4.0 A wall (margin %.3f A)"
           % (name, record["wall"]["initial_margin_angstrom"]),
           record["wall"]["initial_margin_angstrom"] > 0.0)
        ok("%s: HEIGHT literal equals Ha value at 627.509 kcal/mol" % name,
           abs(record["metad"]["height_kcal_mol"]
               - record["metad"]["height_ha"] * prep.HA_TO_KCAL) < 1.0e-7)

    ok("bulk box is 14.2296^3 A", cases["bulk_r1"]["box_angstrom"] == [14.2296] * 3)
    ok("interface box is 15 x 15 x 30 A",
       cases["interface_r1"]["box_angstrom"] == [15.0, 15.0, 30.0])
    iface = cases["interface_r1"]["interface_geometry"]
    ok("interface snapshot is a slab with a %.1f A vacuum gap"
       % iface["vacuum_gap_angstrom"], iface["vacuum_gap_angstrom"] > 10.0)
    ok("interface carbon sits %.2f A below the film top surface"
       % iface["carbon_depth_below_top_surface_angstrom"],
       0.0 < iface["carbon_depth_below_top_surface_angstrom"] < 5.0)
    return manifest, cases


def cp2k_section_balance(text):
    stack = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped.startswith("&"):
            continue
        parts = stripped.split()
        keyword = parts[0][1:].upper()
        if keyword == "END":
            name = parts[1].upper() if len(parts) > 1 else None
            if not stack or stack[-1] != name:
                return "unbalanced &END %s" % name
            stack.pop()
        else:
            stack.append(keyword)
    return None if not stack else "unclosed %s" % stack[-1]


def check_cp2k_inputs(cases):
    for name, record in cases.items():
        text = (HERE / record["inp"]).read_text(encoding="utf-8")
        balance = cp2k_section_balance(text)
        ok("%s: CP2K sections balance" % name, balance is None, str(balance))
        required = [
            "PROJECT " + name,
            "RUN_TYPE MD",
            "CHARGE -2",
            "MULTIPLICITY 1",
            "CUTOFF 300",
            "REL_CUTOFF 40",
            "EPS_SCF 1.0E-6",
            "MINIMIZER DIIS",
            "&XC_FUNCTIONAL BLYP",
            "TYPE DFTD3(BJ)",
            "REFERENCE_FUNCTIONAL BLYP",
            "BASIS_SET DZVP-MOLOPT-SR-GTH",
            "POTENTIAL GTH-BLYP-q1",
            "POTENTIAL GTH-BLYP-q4",
            "ENSEMBLE NVT",
            "TYPE CSVR",
            "TIMECON 100",
            "TIMESTEP 1.0",
            "TEMPERATURE 298.15",
            "STEPS 10000",
            "COORD_FILE_NAME %s.xyz" % name,
            "USE_PLUMED .TRUE.",
            "PLUMED_INPUT_FILE ./%s.plumed" % name,
        ]
        missing = [token for token in required if token not in text]
        ok("%s: CP2K input has all required settings" % name, not missing, str(missing))

        ok("%s: dispersion uses the BJ-damped keyword, not bare DFTD3" % name,
           re.search(r"^\s*TYPE DFTD3\(BJ\)\s*$", text, re.MULTILINE) is not None
           and re.search(r"^\s*TYPE DFTD3\s*$", text, re.MULTILINE) is None)
        ok("%s: CP2K input uses BLYP+D3(BJ) and not PBE" % name,
           "REFERENCE_FUNCTIONAL PBE" not in text and "&XC_FUNCTIONAL PBE" not in text)
        ok("%s: forward GTH-BLYP-q6 appears for both O and S" % name,
           text.count("POTENTIAL GTH-BLYP-q6") == 2)
        box = record["box_angstrom"]
        ok("%s: CELL ABC matches the manifest box" % name,
           "ABC %s %s %s" % (box[0], box[1], box[2]) in text)
        ok("%s: trajectory every 10 steps and restart every 100" % name,
           "FORMAT XYZ" in text and "MD 10" in text and "MD 100" in text)


def check_plumed(cases):
    for name, record in cases.items():
        text = (HERE / record["plumed"]).read_text(encoding="utf-8")
        ok("%s: PLUMED units are A / kcal/mol / fs" % name,
           "UNITS LENGTH=A ENERGY=kcal/mol TIME=fs" in text)
        ok("%s: one CV, d(S-C) with 1-based atoms 293 and 289" % name,
           len(re.findall(r"DISTANCE", text)) == 1
           and "DISTANCE ATOMS=293,289" in text
           and not re.search(r"DISTANCE\s+ATOMS=293,\s+289", text))
        ok("%s: SIGMA 0.1 A and PACE 30 fs" % name,
           "SIGMA=0.1" in text and "PACE=30" in text)
        ok("%s: HEIGHT is the requested kcal/mol literal" % name,
           ("HEIGHT=%s" % record["metad"]["height_kcal_mol"]) in text)
        ok("%s: ordinary metaD has no BIASFACTOR / WT block" % name,
           "BIASFACTOR" not in text.upper() and "WELLTEMPERED" not in text.upper())
        ok("%s: UPPER_WALLS 4.0 A / 200 kcal/mol/A^2" % name,
           "UPPER_WALLS" in text and "AT=4.0" in text and "KAPPA=200.0" in text)
        ok("%s: PRINT d and bias every 10 fs" % name,
           re.search(r"PRINT\s+ARG=d_sc,metad\.bias\s+STRIDE=10", text) is not None)


def check_wrappers():
    wrappers = {
        "bulk_r1": HERE / "submit_bulk_r1.slurm",
        "interface_r1": HERE / "submit_interface_r1.slurm",
    }
    workdirs = set()
    for name, path in wrappers.items():
        text = path.read_text(encoding="utf-8")
        ok("%s: wrapper exists and passes bash -n" % name,
           path.is_file() and subprocess.call(["bash", "-n", str(path)]) == 0)
        ok("%s: 64 MPI ranks, 32G, 48 h" % name,
           "#SBATCH --ntasks=64" in text
           and "#SBATCH --mem=32G" in text
           and "#SBATCH --time=48:00:00" in text)
        ok("%s: compute partition and mail END,FAIL" % name,
           "#SBATCH --partition=compute" in text
           and "#SBATCH --mail-type=END,FAIL" in text)
        ok("%s: no GPU request and one thread per rank" % name,
           "--gres" not in text and "--cpus-per-task" not in text
           and "OMP_NUM_THREADS=1" in text)
        ok("%s: wrapper calls the runner with the matching case" % name,
           'run_cp2k.py" --case %s' % name in text)
        match = re.search(r"^WORK=\$BASE/(\S+)$", text, re.MULTILINE)
        ok("%s: independent fix1 workdir" % name,
           match is not None and match.group(1) == name + "_fix1")
        if match:
            workdirs.add(match.group(1))
    ok("the two wrappers use distinct fix1 workdirs",
       workdirs == {"bulk_r1_fix1", "interface_r1_fix1"})


def check_runner(cases):
    environment = runner.cp2k_environment()
    ok("runner sets OMP_NUM_THREADS=1", environment["OMP_NUM_THREADS"] == "1")
    ok("runner puts the PLUMED shared library on LD_LIBRARY_PATH",
       runner.PLUMED_LIB in environment["LD_LIBRARY_PATH"].split(":"))
    ok("runner points at the PLUMED-enabled cp2k.psmp",
       runner.CP2K.endswith("cp2k.psmp") and "plumed" in runner.CP2K)

    sample = "\n".join([
        " MD| Step number                                           1",
        " SCF run converged in     7 steps",
        " ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]:        -123.456789",
        " MD| Step number                                           2",
        " SCF run converged in     9 steps",
        " ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]:        -123.450000",
        " PROGRAM ENDED AT",
    ])
    parsed = runner.parse_output(sample)
    ok("runner parser reads steps, SCF statistics and energies",
       parsed["program_ended"] and parsed["md_step_first"] == 1
       and parsed["md_step_last"] == 2 and parsed["scf"]["converged_runs"] == 2
       and parsed["scf"]["iterations_max"] == 9
       and parsed["scf"]["iterations_mean"] == 8.0
       and parsed["energies_hartree"]["last"] == -123.45)

    scratch = tempfile.mkdtemp(prefix="li2025_direct_selftest_")
    try:
        trajectory = Path(scratch) / "frames.xyz"
        frame = "3\nframe\nH 0 0 0\nH 0 0 1\nO 0 0 2\n"
        trajectory.write_text(frame + frame, encoding="utf-8")
        ok("runner counts trajectory frames", runner.count_frames(trajectory) == 2)

        for name in ("bulk_r1", "interface_r1"):
            workdir = Path(scratch) / name
            output = subprocess.Popen(
                [sys.executable, str(HERE / "run_cp2k.py"), "--case", name,
                 "--workdir", str(workdir), "--ntasks", "64", "--dry-run",
                 "--workdir-root", scratch],
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0]
            payload = json.loads(output.decode("utf-8"))
            ok("%s: runner dry-run stages inputs and verifies hashes" % name,
               payload["dry_run"] and payload["command"][1:3] == ["-np", "64"]
               and payload["command"][3] == runner.CP2K
               and all((workdir / (name + suffix)).is_file()
                       for suffix in (".xyz", ".inp", ".plumed")))
            ok("%s: runner reports the frozen S-C distance" % name,
               abs(payload["sc_distance_angstrom"]
                   - cases[name]["sc_distance_angstrom"]) < 1.0e-9)

        tampered = Path(scratch) / "bulk_r1" / "bulk_r1.xyz"
        tampered.write_text(tampered.read_text(encoding="utf-8") + "\n", encoding="utf-8")
        code = subprocess.call(
            [sys.executable, str(HERE / "run_cp2k.py"), "--case", "bulk_r1",
             "--workdir", str(Path(scratch) / "bulk_r1"), "--dry-run",
             "--workdir-root", scratch],
            stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        ok("runner refuses a tampered input (exit 4)", code == 4)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def check_python36_and_shell():
    py_files = sorted(HERE.glob("*.py"))
    names = set(path.name for path in py_files)
    required = {"prepare_li2025_direct.py", "run_cp2k.py", "selftest_li2025_direct.py"}
    ok("the first-segment entry points are present", required.issubset(names),
       str(sorted(names)))
    for path in py_files:
        source = path.read_text(encoding="utf-8")
        try:
            ast.parse(source, filename=str(path), feature_version=(3, 6))
            compile(source, str(path), "exec")
        except SyntaxError as error:
            raise AssertionError("%s is not Python 3.6 syntax: %s" % (path.name, error))
        lowered = source
        capture_keyword = "capture" + "_output"
        text_keyword = "te" + "xt"
        for token in (capture_keyword, "missing" + "_ok",
                      "remove" + "prefix", "remove" + "suffix"):
            ok("%s has no %s" % (path.name, token), token not in lowered)
        subprocess_pattern = (r"subprocess\.(run|Popen|call|check_call|check_output)"
                              r"\([^)]*?\b(" + text_keyword + "|" + capture_keyword + r")\s*=")
        ok("%s has no subprocess text/capture keywords" % path.name,
           re.search(subprocess_pattern, source, re.DOTALL) is None)
        ok("%s has no walrus operator" % path.name, (":" + "=") not in source)

    for path in sorted(HERE.glob("*.slurm")):
        ok("%s passes bash -n" % path.name, subprocess.call(["bash", "-n", str(path)]) == 0)


def check_scope_files():
    readme = HERE / "README.md"
    ok("README exists", readme.is_file())
    text = readme.read_text(encoding="utf-8")
    lowered = text.lower()
    for token in ("not li's three-replica", "five-distance", "no fes", "q = -2",
                  "counterion", "project assumption", "controlled"):
        ok("README states %r" % token, token in lowered)
    leftovers = [path.name for path in HERE.rglob("*")
                 if path.suffix.lower() in (".hills", ".colvar")
                 or path.name.startswith(("HILLS_", "COLVAR_"))]
    ok("no HILLS/COLVAR artefacts are present yet (nothing has run)", not leftovers, str(leftovers))


def main():
    exported = check_source_and_export()
    manifest, cases = check_manifest(exported)
    check_cp2k_inputs(cases)
    check_plumed(cases)
    check_wrappers()
    check_runner(cases)
    check_python36_and_shell()
    check_scope_files()
    print("SELFTEST_OK: %d checks passed (source, geometry, CP2K, PLUMED, wrappers, runner, py36, scope)"
          % CHECKS["count"])
    return 0


if __name__ == "__main__":
    sys.exit(main())

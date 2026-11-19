#!/usr/bin/env python3
"""Extract the Li 2025 interface short-contact frame and build three unbiased replicas.

Source (already pulled and SHA-verified):
  results/li2025_direct/interface_r1_seg2/interface_r1-pos-1.xyz   (steps 10010-20000)
  results/li2025_direct/interface_r1_fix1/interface_r1-pos-1.xyz   (steps 0-10000)

The 20 ps interface metaD run (both segments together, 2001 frames) has exactly
six frames with d(S-C) < 2.0 A, all inside one contiguous window:
steps 10370-10420, i.e. 10370-10420 fs. The shortest is step 10410.

Frame identification does NOT compare formatted time strings. The XYZ comment
prints `time = 10410.000` (three decimals) while the COLVAR row for the same step
carries `10409.999160`, so the frame is matched on the MD **step** plus the
numeric **time** inside TIME_TOLERANCE_FS. Zero or several matches stop the build.

What this script does:
  * locates that single frame, copies the 296 atom lines verbatim,
  * verifies composition H194 O100 C1 S1, the fixed atom order, intact waters /
    HCHO / SO3 and the contact gate, all with the 15 x 15 x 30 A interface cell,
  * records the source SHAs, the extracted-frame SHA and the full 20 ps event scan,
  * cross-checks the frame's d(S-C) against the COLVAR row at that step,
  * derives the unbiased inputs from inputs/interface_r1.inp:
        - removes &FREE_ENERGY / &METADYN / USE_PLUMED / PLUMED_INPUT_FILE,
        - adds &GLOBAL SEED=<seed> and no &EXT_RESTART / &VELOCITY,
        - STEPS 3000 (3 ps), COORD_FILE_NAME = the event XYZ,
        - the CELL line must stay `ABC 15.0 15.0 30.0`; the bulk 14.2296 A cube
          is rejected outright,
        - every other scientific line is byte-identical,
  * emits three replica inputs that differ ONLY in the SEED value.

Velocity policy: the event XYZ carries positions only. With no &EXT_RESTART and
no &VELOCITY, CP2K draws fresh Maxwell-Boltzmann velocities at 298.15 K from
&GLOBAL SEED. No metaD velocity is inherited.

Outputs (inside campaign/li2025_direct/):
  inputs/interface_r1_event10410.xyz
  inputs/interface_event_unbiased_s{1,2,3}.inp
  interface_event_unbiased_manifest.json

Run: python3 prepare_interface_event_unbiased.py
"""
import json
import re
from pathlib import Path

import prepare_event_unbiased as bulk

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"

SEG2 = ROOT / "results" / "li2025_direct" / "interface_r1_seg2"
SEG1 = ROOT / "results" / "li2025_direct" / "interface_r1_fix1"
SOURCE = SEG2 / "interface_r1-pos-1.xyz"
SOURCE_SEG1 = SEG1 / "interface_r1-pos-1.xyz"
COLVAR = SEG2 / "COLVAR_interface_r1"
REFERENCE_INPUT = INPUTS / "interface_r1.inp"

EVENT_STEP = 10410
EVENT_TIME_FS = 10409.999160
TIME_TOLERANCE_FS = 0.01
ATOM_COUNT = bulk.ATOM_COUNT
EXPECTED_ORDER = bulk.EXPECTED_ORDER
BOX = (15.0, 15.0, 30.0)
CELL_LINE = "ABC 15.0 15.0 30.0"
EVENT_THRESHOLD = 2.0
EXPECTED_SC = 1.808134
EXPECTED_WINDOW_STEPS = (10370, 10420)
EXPECTED_WINDOW_FRAMES = 6
FRAME_STEM = "interface_r1_event10410"
PROJECT = "interface_r1_unbiased"
STEPS = 3000
TRAJECTORY_EVERY_FS = 10
SEEDS = [("interface_event_unbiased_s1", 20261112),
         ("interface_event_unbiased_s2", 20261113),
         ("interface_event_unbiased_s3", 20261114)]


def scan_sc_series(path, box, threshold=None):
    """Yield (step, time_fs, d(S-C)) for every frame, minimum image in `box`."""
    with open(str(path), "r") as handle:
        while True:
            count_line = handle.readline()
            if not count_line:
                return
            if not count_line.strip():
                continue
            atom_count = int(count_line.strip())
            if atom_count != ATOM_COUNT:
                raise SystemExit("%s has %d atoms, expected %d"
                                 % (path, atom_count, ATOM_COUNT))
            comment = handle.readline()
            carbon = None
            sulfur = None
            for index in range(atom_count):
                parts = handle.readline().split()
                if len(parts) != 4:
                    raise SystemExit("bad atom line in %s: %r" % (path, parts))
                if index == 288:
                    carbon = [float(parts[1]), float(parts[2]), float(parts[3])]
                elif index == 292:
                    sulfur = [float(parts[1]), float(parts[2]), float(parts[3])]
            match = re.search(r"i\s*=\s*(\d+)[,\s]+time\s*=\s*([-0-9.]+)", comment)
            if not match:
                raise SystemExit("frame comment lacks i/time in %s: %r" % (path, comment))
            if carbon is None or sulfur is None:
                raise SystemExit("frame without C or S atoms in %s" % path)
            yield (int(match.group(1)), float(match.group(2)),
                   bulk.distance(sulfur, carbon, box))


def twenty_ps_scan(paths, box, threshold):
    """Every frame of both segments; the short-contact rows and their window."""
    rows = []
    files = []
    for path in paths:
        frames = 0
        first = None
        last = None
        for step, time_fs, sc in scan_sc_series(path, box):
            frames += 1
            if first is None:
                first = step
            last = step
            rows.append({"step": step, "time_fs": time_fs, "sc": sc})
        files.append({
            "path": str(path.relative_to(ROOT)),
            "sha256": bulk.sha256(path),
            "bytes": path.stat().st_size,
            "frames": frames,
            "step_range": [first, last],
        })
    if not rows:
        raise SystemExit("no frames scanned from %s" % paths)
    short = [row for row in rows if row["sc"] < threshold]
    if not short:
        raise SystemExit("no frame below %.2f A in the %.1f ps scan"
                         % (threshold, len(rows) / 10.0))
    steps = [row["step"] for row in short]
    contiguous = all((steps[index + 1] - steps[index])
                     == (steps[1] - steps[0]) for index in range(len(steps) - 1))
    shortest = min(rows, key=lambda row: row["sc"])
    return {
        "threshold_angstrom": threshold,
        "trajectories": files,
        "box_angstrom": list(box),
        "frames_scanned": len(rows),
        "time_span_fs": [rows[0]["time_fs"], rows[-1]["time_fs"]],
        "frames_below_threshold": len(short),
        "below_threshold_rows": short,
        "window_steps": [steps[0], steps[-1]],
        "window_time_fs": [short[0]["time_fs"], short[-1]["time_fs"]],
        "window_is_contiguous": contiguous,
        "window_step_spacing_fs": steps[1] - steps[0] if len(steps) > 1 else None,
        "minimum": shortest,
        "single_window": (contiguous and len(short) == (steps[-1] - steps[0])
                          // (steps[1] - steps[0]) + 1),
    }


def colvar_row(time_fs, tolerance_fs):
    """The COLVAR row whose time matches `time_fs`; d_sc is cross-checked."""
    matches = []
    if not COLVAR.is_file():
        return {"path": None, "note": "COLVAR not found"}
    with open(str(COLVAR), "r") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) < 2:
                continue
            try:
                value = float(parts[0])
                d_sc = float(parts[1])
            except ValueError:
                continue
            if abs(value - time_fs) <= tolerance_fs:
                matches.append((value, d_sc))
    if len(matches) != 1:
        return {"path": str(COLVAR.relative_to(ROOT)),
                "matches": len(matches),
                "note": "expected exactly one COLVAR row at this time"}
    return {
        "path": str(COLVAR.relative_to(ROOT)),
        "sha256": bulk.sha256(COLVAR),
        "time_fs": matches[0][0],
        "d_sc_angstrom": matches[0][1],
        "matches": 1,
    }


def main():
    INPUTS.mkdir(parents=True, exist_ok=True)

    frame = bulk.find_event_frame_by_step(SOURCE, EVENT_STEP, EVENT_TIME_FS,
                                          TIME_TOLERANCE_FS)
    elements, coords = bulk.parse_atoms(frame["atom_lines"])
    geometry = bulk.verify_structure(elements, coords, BOX)
    if abs(geometry["sc_distance_angstrom"] - EXPECTED_SC) > 5.0e-4:
        raise SystemExit("S-C %.6f A does not match the expected %.6f A"
                         % (geometry["sc_distance_angstrom"], EXPECTED_SC))

    scan = twenty_ps_scan([SOURCE_SEG1, SOURCE], BOX, EVENT_THRESHOLD)
    if scan["minimum"]["step"] != EVENT_STEP:
        raise SystemExit("the 20 ps minimum is step %d, not %d"
                         % (scan["minimum"]["step"], EVENT_STEP))
    if abs(scan["minimum"]["sc"] - geometry["sc_distance_angstrom"]) > 1.0e-6:
        raise SystemExit("scan and extracted frame disagree on d(S-C)")
    if scan["frames_below_threshold"] != EXPECTED_WINDOW_FRAMES:
        raise SystemExit("%d frames below %.2f A, expected %d"
                         % (scan["frames_below_threshold"], EVENT_THRESHOLD,
                            EXPECTED_WINDOW_FRAMES))
    if tuple(scan["window_steps"]) != EXPECTED_WINDOW_STEPS:
        raise SystemExit("short-contact window %s != %s"
                         % (scan["window_steps"], list(EXPECTED_WINDOW_STEPS)))
    if not scan["single_window"]:
        raise SystemExit("the short-contact frames are not one contiguous window")

    colvar = colvar_row(EVENT_TIME_FS, TIME_TOLERANCE_FS)
    if colvar.get("matches") != 1:
        raise SystemExit("COLVAR cross-check failed: %s" % colvar)
    if abs(colvar["d_sc_angstrom"] - geometry["sc_distance_angstrom"]) > 5.0e-4:
        raise SystemExit("COLVAR d_sc %.6f A != frame d(S-C) %.6f A"
                         % (colvar["d_sc_angstrom"], geometry["sc_distance_angstrom"]))

    comment = ("Li2025 Li-interface metaD short-contact frame i=%d time=%.6f fs "
               "(the XYZ comment of the source frame prints %.3f; the six-decimal "
               "COLVAR time is used); S-C=%.6f A in the 15x15x30 A interface cell; "
               "verbatim atom lines from %s; q=-2 no counterion; positions only"
               % (frame["step"], EVENT_TIME_FS, frame["time_fs_in_file"],
                  geometry["sc_distance_angstrom"], SOURCE.name))
    xyz_text = "\n".join([str(ATOM_COUNT), comment] + frame["atom_lines"]) + "\n"
    xyz_path = INPUTS / (FRAME_STEM + ".xyz")
    xyz_path.write_text(xyz_text, encoding="utf-8")

    reference_text = REFERENCE_INPUT.read_text(encoding="utf-8")
    if reference_text.count(CELL_LINE) != 1:
        raise SystemExit("reference input must carry exactly one %r cell line" % CELL_LINE)

    cases = []
    for name, seed in SEEDS:
        text, removed = bulk.build_unbiased_input(
            reference_text, seed, FRAME_STEM, project=PROJECT, steps=STEPS,
            reference_project="interface_r1", reference_coord="interface_r1.xyz",
            expected_cell=CELL_LINE)
        if removed != bulk.strip_free_energy(reference_text)[1]:
            raise SystemExit("the removed &FREE_ENERGY block changed")
        inp_path = INPUTS / (name + ".inp")
        inp_path.write_text(text, encoding="utf-8")
        cases.append({
            "name": name,
            "seed": seed,
            "unbiased": True,
            "box_class": "interface",
            "cp2k_project": PROJECT,
            "box_angstrom": list(BOX),
            "sc_distance_angstrom": geometry["sc_distance_angstrom"],
            "stage": [
                {"role": "new", "source": xyz_path.name, "target": xyz_path.name,
                 "sha256": bulk.sha256(xyz_path), "bytes": xyz_path.stat().st_size},
                {"role": "new", "source": inp_path.name, "target": inp_path.name,
                 "sha256": bulk.sha256(inp_path), "bytes": inp_path.stat().st_size},
            ],
            "outputs": {
                "output": "%s.out" % name,
                "trajectory": "%s-pos-1.xyz" % PROJECT,
                "energy": "%s-1.ener" % PROJECT,
                "restart": "%s-1.restart" % PROJECT,
            },
            "steps": STEPS,
            "expected_frames": STEPS // TRAJECTORY_EVERY_FS + 1,
            "fes": {
                "status": "not_applicable",
                "note": "unbiased BOMD; no metadynamics, no FES, no verdict.",
            },
        })

    texts = [(INPUTS / (name + ".inp")).read_text(encoding="utf-8") for name, _ in SEEDS]
    normalised = [re.sub(r"^  SEED \d+\n", "  SEED <seed>\n", text, flags=re.MULTILINE)
                  for text in texts]
    if len(set(normalised)) != 1:
        raise SystemExit("replica inputs differ beyond the SEED line")
    seeds = [seed for _, seed in SEEDS]
    if len(set(seeds)) != len(seeds):
        raise SystemExit("replica seeds are not unique")

    manifest = {
        "schema": "li2025_direct_interface_event_unbiased/v1",
        "session": "li2025-direct-interface-unbiased",
        "purpose": ("Three fully unbiased BOMD replicas started from the single short-contact "
                    "frame of the Li 2025 interface metaD run (step %d, d(S-C) = %.6f A in the "
                    "15 x 15 x 30 A cell), harvested from the biased interface_r1 trajectory."
                    % (EVENT_STEP, geometry["sc_distance_angstrom"])),
        "scope": {
            "what_this_is": [
                "Three independent unbiased NVT BOMD replicas (3 ps each) from one frozen "
                "interface event frame.",
                "CP2K-generated Maxwell-Boltzmann velocities at 298.15 K, seeded per replica.",
                "A stability/relaxation check of a short contact after removing all PLUMED bias.",
            ],
            "what_this_is_not": [
                "NOT spontaneous formation: the start frame came from a biased metaD trajectory.",
                "NOT evidence that the contact is stable; the replicas test whether it survives "
                "without the bias.",
                "NOT a free energy, PMF, barrier or rate.",
                "NOT Li's three-replica statistics and not comparable to the five-distance CV.",
                "NOT a geometry optimisation; the frame is used verbatim.",
            ],
            "structure_conclusion_status": ("Distance facts only; any C-S bonding claim must "
                                            "pass the three-layer gate in AGENTS.md section 7."),
        },
        "source": {
            "trajectory": str(SOURCE.relative_to(ROOT)),
            "sha256": bulk.sha256(SOURCE),
            "bytes": SOURCE.stat().st_size,
            "frame": {
                "step": frame["step"],
                "time_fs": EVENT_TIME_FS,
                "time_fs_printed_in_trajectory": frame["time_fs_in_file"],
                "time_tolerance_fs": TIME_TOLERANCE_FS,
                "match_rule": ("MD step equals %d AND |time - %.6f| <= %.3f fs; formatted "
                               "time strings are never compared, because the XYZ comment "
                               "prints three decimals while COLVAR prints six"
                               % (EVENT_STEP, EVENT_TIME_FS, TIME_TOLERANCE_FS)),
                "comment": frame["comment"],
                "atom_count": frame["atom_count"],
                "line_range_1based": [frame["count_index"] + 1,
                                      frame["comment_index"] + 1 + ATOM_COUNT],
                "unique_match": True,
            },
            "xyz": {
                "path": "inputs/%s.xyz" % FRAME_STEM,
                "sha256": bulk.sha256(xyz_path),
                "bytes": xyz_path.stat().st_size,
            },
            "extraction": ("atom lines copied verbatim from the trajectory; only a new "
                           "comment line was written. No wrapping, rotation, translation "
                           "or bond editing."),
        },
        "colvar_cross_check": colvar,
        "twenty_ps_event_scan": scan,
        "box_angstrom": list(BOX),
        "box_class": "interface",
        "geometry": geometry,
        "atom_order": ("1-288 = 96 waters in H,H,O order; 289 C, 290 O, 291 H, 292 H (HCHO); "
                       "293 S, 294-296 O (SO3)"),
        "solute_indices_1based": bulk.SOLUTE,
        "velocity_strategy": {
            "source": "CP2K Maxwell-Boltzmann velocities at TEMPERATURE 298.15 K",
            "reused_metad_velocities": False,
            "how": ("the event XYZ carries positions only; the input has no &EXT_RESTART "
                    "and no &VELOCITY, so CP2K draws fresh velocities from &GLOBAL SEED"),
            "seeds": [{"case": name, "seed": seed} for name, seed in SEEDS],
        },
        "production": {
            "engine": "CP2K 9.1 (cp2k.psmp, PLUMED build but no PLUMED input)",
            "binary": "/public/home/temp/cp2k-9.1/exe/Linux-x86-64-gfortran-mkl-plumed/cp2k.psmp",
            "functional": "BLYP",
            "dispersion": "DFTD3(BJ), REFERENCE_FUNCTIONAL BLYP",
            "basis": "DZVP-MOLOPT-SR-GTH / GTH-BLYP-q1,q6,q4,q6",
            "cutoff_ry": 300,
            "rel_cutoff_ry": 40,
            "charge": -2,
            "multiplicity": 1,
            "counterion": None,
            "ensemble": "NVT",
            "thermostat": "CSVR",
            "temperature_k": 298.15,
            "timecon_fs": 100,
            "timestep_fs": 1.0,
            "steps": STEPS,
            "steps_fs": STEPS,
            "trajectory_every_fs": TRAJECTORY_EVERY_FS,
            "restart_every_fs": 100,
            "cell": CELL_LINE,
            "mpi_ranks": 64,
            "memory": "32G",
            "max_walltime_hours": 48,
            "partition": "compute (a gpu-partition wrapper is provided for CPU slots)",
            "removed_from_reference": "&FREE_ENERGY / &METADYN / USE_PLUMED / PLUMED_INPUT_FILE",
            "unchanged_from_reference": ("BLYP, DFTD3(BJ)/BLYP, DZVP-MOLOPT-SR-GTH, GTH-BLYP "
                                         "pseudopotentials, CUTOFF 300 / REL_CUTOFF 40, "
                                         "EPS_SCF 1.0E-6, CHARGE -2, MULTIPLICITY 1, "
                                         "ABC 15.0 15.0 30.0, CSVR 298.15 K / TIMECON 100 fs, "
                                         "TIMESTEP 1.0, trajectory every 10 fs, restart every "
                                         "100 fs"),
            "assumptions": [
                "3 ps is a stability check, not a converged sampling window.",
                "The replicas are independent and may be submitted as slots become available.",
                "No counterion is present (q = -2 model), so nothing here is the project's "
                "neutral +2 K+ model.",
            ],
        },
        "verdict": None,
        "cases": cases,
    }
    (HERE / "interface_event_unbiased_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(json.dumps({
        "source_sha256": manifest["source"]["sha256"],
        "frame_step": frame["step"],
        "frame_time_fs": EVENT_TIME_FS,
        "frame_time_printed_in_trajectory": frame["time_fs_in_file"],
        "frame_sha256": manifest["source"]["xyz"]["sha256"],
        "sc_distance_angstrom": geometry["sc_distance_angstrom"],
        "box_angstrom": list(BOX),
        "frames_below_2a_in_20ps": scan["frames_below_threshold"],
        "window_steps": scan["window_steps"],
        "cases": [{"name": name, "seed": seed} for name, seed in SEEDS],
    }, indent=2))


if __name__ == "__main__":
    main()

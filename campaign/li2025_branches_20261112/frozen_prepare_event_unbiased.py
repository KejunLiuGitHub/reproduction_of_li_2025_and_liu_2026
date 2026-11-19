#!/usr/bin/env python3
"""Extract the Li 2025 bulk event frame and build three fully unbiased BOMD inputs.

Source (already pulled and verified):
  results/li2025_direct/bulk_r1_seg2_snapshot_20261108T2101/bulk_r1-pos-1.xyz
  frame i = 12270, time = 12270.000 fs, S-C = 1.778073 A, 296 atoms.

What this script does:
  * requires EXACTLY ONE frame whose time is 12270.000 fs (zero or several -> stop),
  * copies the atom lines verbatim (no wrapping, rotation, translation or bond edit),
  * verifies composition H194 O100 C1 S1 and the fixed order
    C289 / O290 / H291 / H292 / S293 / O294-296,
  * records the source SHA, the extracted-frame SHA, contacts and the event window,
  * derives the unbiased CP2K input from inputs/bulk_r1.inp:
        - removes &FREE_ENERGY/&METADYN/USE_PLUMED/PLUMED_INPUT_FILE entirely,
        - adds &GLOBAL SEED=<seed> and no &EXT_RESTART,
        - STEPS 3000 (3 ps), COORD_FILE_NAME = the event XYZ,
        - every other scientific line is byte-identical,
  * emits three replica inputs that differ ONLY in the SEED value.

Velocity policy: the event XYZ carries positions only. There is no &EXT_RESTART
and no &VELOCITY, so CP2K draws fresh Maxwell-Boltzmann velocities at the input
TEMPERATURE (298.15 K) using &GLOBAL SEED. No metaD velocity is inherited.

Outputs (inside campaign/li2025_direct/):
  inputs/bulk_r1_event12270.xyz
  inputs/event_unbiased_s{1,2,3}.inp
  event_unbiased_manifest.json

The helpers below are reused by `prepare_interface_event_unbiased.py`, which
builds the 15 x 15 x 30 A interface replicas; the cell is a parameter there, so
no bulk constant can leak into an interface build.

Run: python3 prepare_event_unbiased.py
"""
import hashlib
import json
import math
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"
SOURCE = ROOT / "results" / "li2025_direct" / "bulk_r1_seg2_snapshot_20261108T2101" / "bulk_r1-pos-1.xyz"
EVENT_ANALYSIS = ROOT / "results" / "li2025_direct" / "event_analysis" / "BULK_EVENT_12249_12299.json"
REFERENCE_INPUT = INPUTS / "bulk_r1.inp"

EVENT_TIME_FS = 12270.0
EVENT_STEP = 12270
ATOM_COUNT = 296
BOX = (14.2296, 14.2296, 14.2296)
EXPECTED_COMPOSITION = {"H": 194, "O": 100, "C": 1, "S": 1}
EXPECTED_ORDER = ["H", "H", "O"] * 96 + ["C", "O", "H", "H", "S", "O", "O", "O"]
SOLUTE = {"C": 289, "O_carbonyl": 290, "H_formyl_1": 291, "H_formyl_2": 292,
          "S": 293, "O_sulfonate_1": 294, "O_sulfonate_2": 295, "O_sulfonate_3": 296}
FRAME_STEM = "bulk_r1_event12270"
PROJECT = "bulk_r1_unbiased"
STEPS = 3000
SEEDS = [("event_unbiased_s1", 20261109),
         ("event_unbiased_s2", 20261110),
         ("event_unbiased_s3", 20261111)]


def sha256(path):
    digest = hashlib.sha256()
    with open(str(path), "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def distance(a, b, box=None):
    d = [a[i] - b[i] for i in range(3)]
    if box is not None:
        d = [d[i] - box[i] * round(d[i] / box[i]) for i in range(3)]
    return math.sqrt(sum(value * value for value in d))


def find_event_frame(path, time_fs):
    """Return the unique frame whose comment time equals time_fs."""
    matches = []
    with open(str(path), "r") as handle:
        lines = handle.read().splitlines()
    for index, line in enumerate(lines):
        if not line.lstrip().startswith("i ="):
            continue
        match = re.search(r"i\s*=\s*(\d+)[,\s]+time\s*=\s*([-0-9.]+)", line)
        if not match:
            continue
        if abs(float(match.group(2)) - time_fs) < 1.0e-9:
            matches.append((index, int(match.group(1)), line))
    if len(matches) == 0:
        raise SystemExit("no frame at time %.3f fs in %s" % (time_fs, path))
    if len(matches) > 1:
        raise SystemExit("%d frames match time %.3f fs in %s; refusing to guess"
                         % (len(matches), time_fs, path))
    comment_index, step, comment = matches[0]
    count_index = comment_index - 1
    if count_index < 0:
        raise SystemExit("frame at time %.3f fs has no atom-count line" % time_fs)
    try:
        atom_count = int(lines[count_index].strip())
    except ValueError:
        raise SystemExit("bad atom-count line before frame at time %.3f fs" % time_fs)
    if atom_count != ATOM_COUNT:
        raise SystemExit("frame at time %.3f fs has %d atoms, expected %d"
                         % (time_fs, atom_count, ATOM_COUNT))
    atom_lines = lines[comment_index + 1:comment_index + 1 + atom_count]
    if len(atom_lines) != atom_count:
        raise SystemExit("frame at time %.3f fs is truncated" % time_fs)
    return {
        "count_index": count_index,
        "comment_index": comment_index,
        "atom_lines": atom_lines,
        "comment": comment,
        "step": step,
        "atom_count": atom_count,
    }


def parse_atoms(atom_lines):
    elements = []
    coords = []
    for offset, line in enumerate(atom_lines):
        parts = line.split()
        if len(parts) != 4 or parts[0] not in ("H", "C", "O", "S"):
            raise SystemExit("unparseable atom line %d: %r" % (offset + 1, line))
        try:
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
        except ValueError:
            raise SystemExit("unparseable coordinate at atom line %d: %r" % (offset + 1, line))
        elements.append(parts[0])
    return elements, coords


def verify_structure(elements, coords, box=BOX):
    """Composition/order/connectivity/contact gate; `box` sets the minimum image.

    The interface event build calls this with box (15, 15, 30); the bulk build
    keeps the 14.2296 A cube.
    """
    if len(elements) != ATOM_COUNT:
        raise SystemExit("expected %d atoms, got %d" % (ATOM_COUNT, len(elements)))
    composition = dict(Counter(elements))
    if composition != EXPECTED_COMPOSITION:
        raise SystemExit("composition %s != %s" % (composition, EXPECTED_COMPOSITION))
    if elements != EXPECTED_ORDER:
        raise SystemExit("atom order differs from 96x(H,H,O) + C,O,H,H,S,O,O,O")

    # Duplicates and severe overlaps, minimum image.
    closest = (1.0e9, None)
    for i in range(ATOM_COUNT):
        for j in range(i + 1, ATOM_COUNT):
            d = distance(coords[i], coords[j], box)
            if d < closest[0]:
                closest = (d, (i + 1, j + 1, elements[i], elements[j]))
    if closest[0] < 1.0e-6:
        raise SystemExit("exact duplicate atoms: %s (%.3e A)" % (closest[1], closest[0]))
    if closest[0] < 0.6:
        raise SystemExit("severe overlap at %.3f A: %s" % (closest[0], closest[1]))

    water_atoms = []
    for start in range(0, 288, 3):
        water_atoms.append([start, start + 1, start + 2])
    for water in water_atoms:
        oxygen = water[2]
        for hydrogen in (water[0], water[1]):
            if distance(coords[oxygen], coords[hydrogen], box) > 1.30:
                raise SystemExit("water at 1-based O%d is not intact" % (oxygen + 1))
    for oxygen in (293, 294, 295):
        if distance(coords[292], coords[oxygen], box) > 1.75:
            raise SystemExit("SO3 is not intact (S-O%d)" % (oxygen + 1))
    if distance(coords[288], coords[289], box) > 1.45:
        raise SystemExit("HCHO C=O is not intact")
    for hydrogen in (290, 291):
        if distance(coords[288], coords[hydrogen], box) > 1.25:
            raise SystemExit("HCHO C-H is not intact")

    groups = [water for water in water_atoms] + [[288, 289, 290, 291], [292, 293, 294, 295]]
    labels = ["H2O"] * 96 + ["HCHO", "SO3"]
    best = (1.0e9, None)
    for a in range(len(groups)):
        for b in range(a + 1, len(groups)):
            for i in groups[a]:
                for j in groups[b]:
                    d = distance(coords[i], coords[j], box)
                    if d < best[0]:
                        best = (d, (i + 1, labels[a], j + 1, labels[b],
                                    elements[i], elements[j]))
    if best[0] < 1.20:
        raise SystemExit("severe intermolecular overlap at %.3f A: %s" % (best[0], best[1]))

    return {
        "composition": composition,
        "min_atom_pair_distance_including_bonds_angstrom": closest[0],
        "min_atom_pair_including_bonds_1based": list(closest[1]),
        "min_intermolecular_contact_angstrom": best[0],
        "min_intermolecular_contact_atoms_1based": list(best[1]),
        "sc_distance_angstrom": distance(coords[292], coords[288], box),
        "co_distance_angstrom": distance(coords[288], coords[289], box),
        "ch_distances_angstrom": [distance(coords[288], coords[290], box),
                                  distance(coords[288], coords[291], box)],
        "so_distances_angstrom": [distance(coords[292], coords[oxygen], box)
                                  for oxygen in (293, 294, 295)],
    }


def strip_free_energy(text):
    lines = text.splitlines(keepends=True)
    starts = [index for index, line in enumerate(lines)
              if line.strip() == "&FREE_ENERGY"]
    ends = [index for index, line in enumerate(lines)
            if line.strip() == "&END FREE_ENERGY"]
    if len(starts) != 1 or len(ends) != 1 or ends[0] <= starts[0]:
        raise SystemExit("reference input must contain exactly one &FREE_ENERGY block")
    block = [line.rstrip("\n") for line in lines[starts[0]:ends[0] + 1]]
    del lines[starts[0]:ends[0] + 1]
    return "".join(lines), block


def replace_exactly_once(text, old, new, label):
    count = text.count(old)
    if count != 1:
        raise SystemExit("expected exactly one %r (%s), found %d" % (old, label, count))
    return text.replace(old, new, 1)


def build_unbiased_input(reference_text, seed, frame_name, project=PROJECT,
                         steps=STEPS, reference_project="bulk_r1",
                         reference_coord="bulk_r1.xyz", expected_cell=None):
    """Derive one unbiased input; every scientific line stays byte-identical.

    The defaults reproduce the bulk event build. The interface event build
    passes its own project/frame name and `expected_cell="ABC 15.0 15.0 30.0"`
    so a 14.2296 A cube can never be emitted by mistake.
    """
    text, removed = strip_free_energy(reference_text)
    text = replace_exactly_once(text, "  PROJECT %s\n" % reference_project,
                                "  PROJECT %s\n" % project, "PROJECT")
    text, count = re.subn(r"^(\s*)STEPS 10000\s*$", r"\g<1>STEPS %d" % steps,
                          text, count=1, flags=re.MULTILINE)
    if count != 1:
        raise SystemExit("expected exactly one STEPS 10000 line, found %d" % count)
    text = replace_exactly_once(text, "      COORD_FILE_NAME %s\n" % reference_coord,
                                "      COORD_FILE_NAME %s.xyz\n" % frame_name, "COORD")
    text = replace_exactly_once(text, "  PRINT_LEVEL MEDIUM\n",
                                "  PRINT_LEVEL MEDIUM\n  SEED %d\n" % seed, "SEED")
    upper = text.upper()
    for forbidden in ("&FREE_ENERGY", "&METADYN", "USE_PLUMED", "PLUMED_INPUT_FILE",
                      "&EXT_RESTART", "PLUMED"):
        if forbidden in upper:
            raise SystemExit("unbiased input still contains %s" % forbidden)
    if text.count("SEED ") != 1:
        raise SystemExit("unbiased input must contain exactly one SEED line")
    if text.count("STEPS %d" % steps) != 1:
        raise SystemExit("unbiased input lacks the single STEPS %d line" % steps)
    if expected_cell is not None:
        if text.count(expected_cell) != 1:
            raise SystemExit("unbiased input lacks the single %r cell line" % expected_cell)
        if "14.2296" in text:
            raise SystemExit("unbiased input still carries the bulk 14.2296 A cell")
    return text, removed


def find_event_frame_by_step(path, step, time_fs, tolerance_fs):
    """Return the unique frame matching BOTH the MD step and the numeric time.

    The XYZ comment prints time with three decimals while the trailing COLVAR /
    report value carries six, so the time is matched numerically inside
    `tolerance_fs` and never by comparing formatted strings. The step makes the
    match unambiguous: neighbouring frames are one MD step apart.
    """
    matches = []
    with open(str(path), "r") as handle:
        lines = handle.read().splitlines()
    for index, line in enumerate(lines):
        if not line.lstrip().startswith("i ="):
            continue
        match = re.search(r"i\s*=\s*(\d+)[,\s]+time\s*=\s*([-0-9.]+)", line)
        if not match:
            continue
        if int(match.group(1)) != step:
            continue
        if abs(float(match.group(2)) - time_fs) > tolerance_fs:
            continue
        matches.append((index, int(match.group(1)), float(match.group(2)), line))
    if len(matches) == 0:
        raise SystemExit("no frame with step %d and time %.6f +/- %.6f fs in %s"
                         % (step, time_fs, tolerance_fs, path))
    if len(matches) > 1:
        raise SystemExit("%d frames match step %d and time %.6f fs in %s; refusing to guess"
                         % (len(matches), step, time_fs, path))
    comment_index, matched_step, matched_time, comment = matches[0]
    count_index = comment_index - 1
    if count_index < 0:
        raise SystemExit("frame at step %d has no atom-count line" % step)
    try:
        atom_count = int(lines[count_index].strip())
    except ValueError:
        raise SystemExit("bad atom-count line before frame at step %d" % step)
    if atom_count != ATOM_COUNT:
        raise SystemExit("frame at step %d has %d atoms, expected %d"
                         % (step, atom_count, ATOM_COUNT))
    atom_lines = lines[comment_index + 1:comment_index + 1 + atom_count]
    if len(atom_lines) != atom_count:
        raise SystemExit("frame at step %d is truncated" % step)
    return {
        "count_index": count_index,
        "comment_index": comment_index,
        "atom_lines": atom_lines,
        "comment": comment,
        "step": matched_step,
        "time_fs_in_file": matched_time,
        "atom_count": atom_count,
    }


def scan_sc_min(path, box=BOX):
    """Global minimum d(S-C) across the whole trajectory, with its step/time."""
    best = (1.0e9, None, None)
    with open(str(path), "r") as handle:
        while True:
            count_line = handle.readline()
            if not count_line:
                break
            if not count_line.strip():
                continue
            atom_count = int(count_line.strip())
            comment = handle.readline()
            carbon = None
            sulfur = None
            for index in range(atom_count):
                parts = handle.readline().split()
                if index == 288:
                    carbon = [float(parts[1]), float(parts[2]), float(parts[3])]
                elif index == 292:
                    sulfur = [float(parts[1]), float(parts[2]), float(parts[3])]
            if carbon is None or sulfur is None:
                continue
            value = distance(sulfur, carbon, box)
            if value < best[0]:
                match = re.search(r"i\s*=\s*(\d+)[,\s]+time\s*=\s*([-0-9.]+)", comment)
                best = (value,
                        int(match.group(1)) if match else None,
                        float(match.group(2)) if match else None)
    return {"sc_angstrom": best[0], "step": best[1], "time_fs": best[2]}


def event_context():
    if not EVENT_ANALYSIS.is_file():
        return {"path": None, "note": "event analysis JSON not found"}
    data = json.loads(EVENT_ANALYSIS.read_text(encoding="utf-8"))
    rows = data.get("rows", [])
    sc_at_event = None
    window_min = None
    for row in rows:
        if abs(float(row.get("time_fs", -1.0)) - EVENT_TIME_FS) < 1.0e-9:
            sc_at_event = row.get("SC")
        value = row.get("SC")
        if value is not None and (window_min is None or value < window_min):
            window_min = value
    return {
        "path": str(EVENT_ANALYSIS.relative_to(ROOT)),
        "sha256": sha256(EVENT_ANALYSIS),
        "window_time_fs": [rows[0]["time_fs"], rows[-1]["time_fs"]] if rows else None,
        "sc_at_event_from_analysis_angstrom": sc_at_event,
        "event_window_min_sc_angstrom": window_min,
        "note": data.get("note"),
    }


def _with_global_sc_context(context):
    global_min = scan_sc_min(SOURCE)
    context["trajectory_min_sc_angstrom"] = global_min["sc_angstrom"]
    context["trajectory_min_sc_step"] = global_min["step"]
    context["trajectory_min_sc_time_fs"] = global_min["time_fs"]
    context["requested_frame_is_global_min"] = (
        global_min["step"] == EVENT_STEP
        and abs(global_min["sc_angstrom"] - 1.778073) < 1.0e-6)
    context["note"] = (
        "The 12270 fs frame is the event-window minimum and the frame requested for this "
        "validation. The trajectory reaches a slightly shorter S-C (%0.4f A at step %s, "
        "%0.1f fs); that frame is NOT used here."
        % (global_min["sc_angstrom"], global_min["step"], global_min["time_fs"]))
    return context


def main():
    INPUTS.mkdir(parents=True, exist_ok=True)

    frame = find_event_frame(SOURCE, EVENT_TIME_FS)
    elements, coords = parse_atoms(frame["atom_lines"])
    geometry = verify_structure(elements, coords)
    if abs(geometry["sc_distance_angstrom"] - 1.778073) > 5.0e-4:
        raise SystemExit("S-C %.6f A does not match the expected 1.778073 A"
                         % geometry["sc_distance_angstrom"])

    comment = ("Li2025 Li-bulk metaD event frame i=%d time=%.3f fs; S-C=%.6f A; "
               "verbatim from %s; q=-2 no K/Na; positions only"
               % (frame["step"], EVENT_TIME_FS, geometry["sc_distance_angstrom"],
                  SOURCE.name))
    xyz_text = "\n".join([str(ATOM_COUNT), comment] + frame["atom_lines"]) + "\n"
    xyz_path = INPUTS / (FRAME_STEM + ".xyz")
    xyz_path.write_text(xyz_text, encoding="utf-8")

    reference_text = REFERENCE_INPUT.read_text(encoding="utf-8")
    cases = []
    for name, seed in SEEDS:
        text, removed = build_unbiased_input(reference_text, seed, FRAME_STEM)
        inp_path = INPUTS / (name + ".inp")
        inp_path.write_text(text, encoding="utf-8")
        cases.append({
            "name": name,
            "seed": seed,
            "unbiased": True,
            "cp2k_project": PROJECT,
            "box_angstrom": list(BOX),
            "sc_distance_angstrom": geometry["sc_distance_angstrom"],
            "stage": [
                {"role": "new", "source": xyz_path.name, "target": xyz_path.name,
                 "sha256": sha256(xyz_path), "bytes": xyz_path.stat().st_size},
                {"role": "new", "source": inp_path.name, "target": inp_path.name,
                 "sha256": sha256(inp_path), "bytes": inp_path.stat().st_size},
            ],
            "outputs": {
                "output": "%s.out" % name,
                "trajectory": "%s-pos-1.xyz" % PROJECT,
                "energy": "%s-1.ener" % PROJECT,
                "restart": "%s-1.restart" % PROJECT,
            },
            "steps": STEPS,
            "expected_frames": STEPS // 10 + 1,
            "fes": {
                "status": "not_applicable",
                "note": "unbiased BOMD; no metadynamics, no FES.",
            },
        })

    # The three inputs must differ only in the SEED value.
    texts = [(HERE / "inputs" / (name + ".inp")).read_text(encoding="utf-8")
             for name, _ in SEEDS]
    normalised = [re.sub(r"^  SEED \d+\n", "  SEED <seed>\n", text, flags=re.MULTILINE)
                  for text in texts]
    if len(set(normalised)) != 1:
        raise SystemExit("replica inputs differ beyond the SEED line")
    seeds = [seed for _, seed in SEEDS]
    if len(set(seeds)) != len(seeds):
        raise SystemExit("replica seeds are not unique")

    source_sha = sha256(SOURCE)
    manifest = {
        "schema": "li2025_direct_event_unbiased/v1",
        "session": "li2025-direct-unbiased-validation",
        "purpose": ("First fully unbiased BOMD replicas started from the requested Li 2025 "
                    "bulk event frame at 12270 fs (S-C = 1.778073 A, the minimum inside "
                    "the 12220-12320 fs event window), harvested from the biased bulk_r1 "
                    "metaD trajectory."),
        "scope": {
            "what_this_is": [
                "Three independent unbiased NVT BOMD replicas (3 ps each) from one frozen event frame.",
                "CP2K-generated Maxwell-Boltzmann velocities at 298.15 K, seeded per replica.",
                "A stability/relaxation check after removing all PLUMED bias.",
            ],
            "what_this_is_not": [
                "NOT spontaneous formation: the starting frame came from a metaD trajectory.",
                "NOT a free energy, PMF, barrier or rate.",
                "NOT Li's three-replica statistics and not comparable to the five-distance CV.",
                "NOT a geometry optimisation; the frame is used verbatim.",
            ],
            "structure_conclusion_status": ("Distance facts only; any C-S bonding claim must "
                                            "pass the three-layer gate in AGENTS.md section 7."),
        },
        "source": {
            "trajectory": str(SOURCE.relative_to(ROOT)),
            "sha256": source_sha,
            "bytes": SOURCE.stat().st_size,
            "frame": {
                "step": frame["step"],
                "time_fs": EVENT_TIME_FS,
                "comment": frame["comment"],
                "atom_count": frame["atom_count"],
                "line_range_1based": [frame["count_index"] + 1,
                                      frame["comment_index"] + 1 + ATOM_COUNT],
                "unique_match": True,
            },
            "xyz": {
                "path": "inputs/%s.xyz" % FRAME_STEM,
                "sha256": sha256(xyz_path),
                "bytes": xyz_path.stat().st_size,
            },
            "extraction": ("atom lines copied verbatim from the trajectory; only a new "
                           "comment line was written. No wrapping, rotation, translation "
                           "or bond editing."),
        },
        "geometry": geometry,
        "atom_order": ("1-288 = 96 waters in H,H,O order; 289 C, 290 O, 291 H, 292 H (HCHO); "
                       "293 S, 294-296 O (SO3)"),
        "solute_indices_1based": SOLUTE,
        "event": _with_global_sc_context(event_context()),
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
            "ensemble": "NVT",
            "thermostat": "CSVR",
            "temperature_k": 298.15,
            "timecon_fs": 100,
            "timestep_fs": 1.0,
            "steps": STEPS,
            "steps_fs": STEPS,
            "trajectory_every_fs": 10,
            "restart_every_fs": 100,
            "mpi_ranks": 64,
            "memory": "32G",
            "max_walltime_hours": 48,
            "partition": "compute",
            "removed_from_reference": "&FREE_ENERGY / &METADYN / USE_PLUMED / PLUMED_INPUT_FILE",
            "assumptions": [
                "3 ps is a first stability check, not a converged sampling window.",
                "CSVR TIMECON 100 fs is unchanged from the production reference.",
                "The three replicas may be submitted as slots become available; they are "
                "independent and need not start together.",
            ],
        },
        "cases": cases,
    }
    (HERE / "event_unbiased_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(json.dumps({
        "source_sha256": source_sha,
        "frame_step": frame["step"],
        "frame_time_fs": EVENT_TIME_FS,
        "frame_sha256": sha256(xyz_path),
        "sc_distance_angstrom": geometry["sc_distance_angstrom"],
        "min_intermolecular_contact_angstrom": geometry["min_intermolecular_contact_angstrom"],
        "cases": [{"name": name, "seed": seed} for name, seed in SEEDS],
    }, indent=2))


if __name__ == "__main__":
    main()

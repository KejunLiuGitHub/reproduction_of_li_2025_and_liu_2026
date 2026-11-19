#!/usr/bin/env python3
"""Shared helpers for the unbiased-interface B3LYP RC seed pipeline.

The interface is a periodic 17 x 17 x 47 A water slab (128 H2O) with neutral
HCHO + HOSO2- + K+ at one free surface.  This module only handles the
*isolated 13-atom RC cluster* that is optimised in the gas phase with SMD
water and then rigidly embedded back into one unbiased MD frame.

Reuse:
  * ``analyze_iface.read_frames`` / ``frame_metrics`` for the strict 394-atom
    trajectory parse and the slab/interface metrics.
  * ``analyze_unbiased.load_trajectory`` for step/time parsing, atom-count and
    element checks, deduplication and monotonically increasing time.

Python 3.6 compatible: the cluster interpreter is old, so no
``subprocess.run(capture_output=...)``, no dataclasses, no walrus.

The 13 cluster atoms keep the source order; nothing is reordered:

  local 0-3   HCHO        source 385-388 (0-based 384-387)
  local 4-8   HOSO2-      source 389-393 (0-based 388-392)
  local 9     K+          source 394   (0-based 393)
  local 10-12 mediator    source water O,H,H (0-based 3w, 3w+1, 3w+2)
"""
import hashlib
import json
import os
import re
import sys
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import analyze_iface  # noqa: E402
import analyze_unbiased  # noqa: E402


EXPECTED_ATOMS = 394
WATER_COUNT = 128
CELL = tuple(float(value) for value in analyze_iface.CELL)

# 0-based source indices.
IDX_C = 384
IDX_O_CARBONYL = 385
IDX_H_ALDEHYDE_1 = 386
IDX_H_ALDEHYDE_2 = 387
IDX_S = 388
IDX_O_HYDROXYL = 389
IDX_O_NONHYDROXYL_1 = 390
IDX_O_NONHYDROXYL_2 = 391
IDX_H_HYDROXYL = 392
IDX_K = 393

SOLUTE_SOURCE_INDICES = (
    IDX_C, IDX_O_CARBONYL, IDX_H_ALDEHYDE_1, IDX_H_ALDEHYDE_2,
    IDX_S, IDX_O_HYDROXYL, IDX_O_NONHYDROXYL_1, IDX_O_NONHYDROXYL_2,
    IDX_H_HYDROXYL, IDX_K,
)
SOLUTE_ROLES = (
    "HCHO_C", "HCHO_O_carbonyl", "HCHO_H1", "HCHO_H2",
    "HOSO2_S", "HOSO2_O_hydroxyl", "HOSO2_O_nonhydroxyl_1",
    "HOSO2_O_nonhydroxyl_2", "HOSO2_H_hydroxyl", "K_counterion",
)
SOLUTE_ELEMENTS = ("C", "O", "H", "H", "S", "O", "O", "O", "H", "K")

# Local indices in the 13-atom cluster.
LOCAL_C = 0
LOCAL_O_CARBONYL = 1
LOCAL_S = 4
LOCAL_O_HYDROXYL = 5
LOCAL_O_NONHYDROXYL_1 = 6
LOCAL_O_NONHYDROXYL_2 = 7
LOCAL_H_HYDROXYL = 8
LOCAL_K = 9
LOCAL_MEDIATOR_O = 10
LOCAL_MEDIATOR_H1 = 11
LOCAL_MEDIATOR_H2 = 12
CLUSTER_SIZE = 13

# All non-hydrogen atoms of the 13-atom cluster (C, O386, S, O390-O392, K, Ow).
HEAVY_LOCAL_INDICES = (0, 1, 4, 5, 6, 7, 9, 10)

LEVEL_LINE = "! B3LYP D3BJ 6-311+G(d,p) Opt NumFreq TightSCF miniprint SMD(water)"
MAXCORE = 4000       # matches campaign/viewer/orca_parser.orca_input()
NPROCS = 1           # ditto; <=20 atoms is fastest on one core
CHARGE = 0
MULTIPLICITY = 1

IMAGINARY_FAIL_BELOW_CM1 = -10.0

# Candidate criteria (project decision, not literature thresholds).
SC_MIN_A = 2.6
SC_MAX_A = 3.8
OH_MAX_A = 1.3
NONHYDROXYL_OC_MIN_A = 2.4
BRIDGE_MAX_A = 2.2
MEDIATOR_OH_MAX_A = 1.3
INTERFACE_MAX_A = 5.0

# Post-optimisation / embedding criteria.
POST_SC_MIN_A = 2.4
POST_OH_MAX_A = 1.3
POST_NONHYDROXYL_OC_MIN_A = 2.4
EMBED_BRIDGE_MAX_A = 2.5
EMBED_INTERFACE_MAX_A = 5.0
EMBED_MIN_EXTERNAL_A = 1.4

REASON_ORDER = (
    "outside_last_5ps",
    "frame_metrics_failed",
    "sc_out_of_window",
    "hydroxyl_oh_not_intact",
    "nonhydroxyl_o_c_too_close",
    "no_mediating_water",
    "not_near_interface",
)

ORCA_INPUT_TEMPLATE = """{level}
%maxcore {maxcore}
%pal nprocs {nprocs} end
%scf MaxIter 800 end
%geom
  MaxIter 300
{constraints}
end
* xyzfile {charge} {mult} {xyz_name}
"""

ENERGY_RE = re.compile(
    r"FINAL SINGLE POINT ENERGY\s+(-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)")
COORD_BLOCK_RE = re.compile(
    r"CARTESIAN COORDINATES \(ANGSTROEM\)\s*\n\s*-+\s*\n(.*?)(?:\n\s*\n|\Z)",
    re.S)
FREQ_LINE_RE = re.compile(
    r"^\s*(\d+):\s+(-?(?:\d+(?:\.\d*)?|\.\d+))\s+cm\*\*-1"
    r"(?:\s+\*\*\*imaginary mode\*\*\*)?\s*$", re.M)
XYZFILE_RE = re.compile(
    r"^\s*\*\s*xyzfile\s+(-?\d+)\s+(\d+)\s+(\S+)\s*$", re.M)
ELEMENT_RE = re.compile(r"^[A-Z][a-z]?$")


class RcError(RuntimeError):
    """Raised when an input or output violates the documented contract."""


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path):
    digest = hashlib.sha256()
    with open(str(path), "rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def sha256_text(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def atomic_write_text(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix="." + path.name + ".", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, str(path))
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def atomic_write_json(path, data):
    atomic_write_text(path, json.dumps(data, indent=2, sort_keys=False) + "\n")


def resolve_against(base, value):
    path = Path(value)
    return path if path.is_absolute() else (Path(base) / path)


def canonical_atoms_text(atoms):
    """Canonical, order-preserving text used for frame hashes."""
    return "".join(
        "{} {:.10f} {:.10f} {:.10f}\n".format(atom[0], atom[1], atom[2], atom[3])
        for atom in atoms)


def canonical_atoms_sha256(atoms):
    return sha256_text(canonical_atoms_text(atoms))


def xyz_text(atoms, comment):
    body = "".join(
        "{:<2s} {:16.10f} {:16.10f} {:16.10f}\n".format(
            atom[0], atom[1], atom[2], atom[3])
        for atom in atoms)
    return "{}\n{}\n{}".format(len(atoms), comment, body)


def write_xyz(path, atoms, comment):
    atomic_write_text(path, xyz_text(atoms, comment))


def read_xyz(path):
    """Read one XYZ frame (strict, no partial block)."""
    lines = Path(path).read_text().splitlines()
    while lines and not lines[-1].strip():
        lines.pop()
    if len(lines) < 2:
        raise RcError("{}: not an XYZ file".format(path))
    try:
        count = int(lines[0].strip())
    except ValueError:
        raise RcError("{}: bad XYZ atom count {!r}".format(path, lines[0]))
    if count <= 0 or len(lines) != count + 2:
        raise RcError("{}: expected {} atom lines, found {}".format(
            path, count, len(lines) - 2))
    atoms = []
    for line_number, line in enumerate(lines[2:], 3):
        fields = line.split()
        if len(fields) < 4:
            raise RcError("{}: line {}: malformed atom line {!r}".format(
                path, line_number, line))
        if not ELEMENT_RE.match(fields[0]):
            raise RcError("{}: line {}: bad element {!r}".format(
                path, line_number, fields[0]))
        try:
            coords = (float(fields[1]), float(fields[2]), float(fields[3]))
        except ValueError:
            raise RcError("{}: line {}: non-numeric coordinate {!r}".format(
                path, line_number, line))
        atoms.append((fields[0],) + coords)
    return atoms, lines[1]


def periodic_distance(first, second, box=CELL):
    return analyze_unbiased.periodic_distance(
        tuple(first), tuple(second), tuple(float(value) for value in box))


def check_source_atoms(atoms, path, frame_start_line):
    if len(atoms) != EXPECTED_ATOMS:
        raise RcError(
            "{}: frame at line {} has {} atoms, expected {}".format(
                path, frame_start_line, len(atoms), EXPECTED_ATOMS))
    for local, source_index in enumerate(SOLUTE_SOURCE_INDICES):
        element = str(atoms[source_index][0]).strip()
        expected = SOLUTE_ELEMENTS[local]
        if element.upper() != expected:
            raise RcError(
                "{}: frame at line {}: source atom {} must be {}, found {!r}".format(
                    path, frame_start_line, source_index + 1, expected, element))
    for water in range(WATER_COUNT):
        base = 3 * water
        found = tuple(str(atoms[base + offset][0]).strip().upper() for offset in range(3))
        if found != ("O", "H", "H"):
            raise RcError(
                "{}: frame at line {}: water {} atoms {} must be O,H,H, found "
                "{}".format(path, frame_start_line, water + 1, (base + 1, base + 2,
                                                                base + 3), found))


def load_prepared(path):
    """Read IFACE_PREPARED.json and validate the fields this pipeline uses."""
    path = Path(path)
    if not path.is_file():
        raise RcError("prepared manifest not found: {}".format(path))
    try:
        data = json.loads(path.read_text())
    except ValueError as exc:
        raise RcError("{}: not valid JSON ({})".format(path, exc))
    if data.get("status") not in (None, "prepared"):
        raise RcError("{}: status is {!r}, expected 'prepared'".format(
            path, data.get("status")))
    if int(data.get("atom_count", -1)) != EXPECTED_ATOMS:
        raise RcError("{}: atom_count must be {}, found {!r}".format(
            path, EXPECTED_ATOMS, data.get("atom_count")))
    cell = data.get("cell_angstrom")
    if not isinstance(cell, (list, tuple)) or len(cell) != 3:
        raise RcError("{}: cell_angstrom must be three lengths".format(path))
    for actual, expected in zip(cell, CELL):
        if abs(float(actual) - expected) > 1.0e-8:
            raise RcError(
                "{}: cell {} does not match the analyze_iface cell {}".format(
                    path, tuple(cell), CELL))
    trajectory_field = data.get("trajectory_path")
    if not trajectory_field:
        trajectory = (data.get("trajectory") or {}).get("path")
        if not trajectory:
            raise RcError("{}: no trajectory_path".format(path))
        trajectory_field = trajectory
    base = path.resolve().parent
    trajectory = resolve_against(base, trajectory_field)
    if not trajectory.is_file():
        raise RcError("source trajectory not found: {}".format(trajectory))
    data["_manifest_path"] = str(path.resolve())
    data["_manifest_sha256"] = sha256_file(path)
    data["_trajectory_path"] = str(trajectory.resolve())
    data["_trajectory_sha256"] = sha256_file(trajectory)
    return data


def load_source_trajectory(path):
    """Strict 394-atom load, then the extra water/element checks.

    Returns (records, duplicate_info, start_lines).  ``start_lines`` maps the
    deduplicated step to its frame start line, with the last occurrence in the
    file winning (the same policy as analyze_unbiased).
    """
    path = Path(path)
    try:
        records, duplicate_info = analyze_unbiased.load_trajectory(path, CELL)
    except Exception as exc:
        raise RcError("{}: strict trajectory parse failed: {}".format(path, exc))
    frames = analyze_iface.read_frames(path)
    start_lines = {}
    for frame in frames:
        if frame.step is not None:
            start_lines[frame.step] = frame.start_line
    for record in records:
        check_source_atoms(record.atoms, path, start_lines[record.step])
    return records, duplicate_info, start_lines


def evaluate_frame(record, start_line, box=CELL):
    """Evaluate one frame against every candidate criterion.

    Never raises for a merely bad frame: the reason is recorded instead.  A
    ``frame_metrics`` failure is one rejection reason, not a crash, so the
    rejection statistics stay complete.
    """
    atoms = record.atoms
    c = atoms[IDX_C][1:4]
    o_carbonyl = atoms[IDX_O_CARBONYL][1:4]
    s = atoms[IDX_S][1:4]
    o_hydroxyl = atoms[IDX_O_HYDROXYL][1:4]
    h_hydroxyl = atoms[IDX_H_HYDROXYL][1:4]
    o_non_1 = atoms[IDX_O_NONHYDROXYL_1][1:4]
    o_non_2 = atoms[IDX_O_NONHYDROXYL_2][1:4]

    sc = periodic_distance(s, c, box)
    oh = periodic_distance(o_hydroxyl, h_hydroxyl, box)
    oc1 = periodic_distance(o_non_1, c, box)
    oc2 = periodic_distance(o_non_2, c, box)
    nonhydroxyl_min = min(oc1, oc2)

    mediator = find_mediator(atoms, box=box)
    metrics = None
    frame_metrics_error = None
    try:
        metrics = analyze_iface.frame_metrics(atoms).result
    except Exception as exc:                      # report, do not hide
        frame_metrics_error = str(exc)

    d_iface_c = None
    d_iface_s = None
    if metrics is not None:
        d_iface_c = metrics["solutes"]["C385"]["distance_to_nearest_interface"]
        d_iface_s = metrics["solutes"]["S389"]["distance_to_nearest_interface"]

    reasons = []
    if frame_metrics_error is not None:
        reasons.append("frame_metrics_failed")
    if not (SC_MIN_A <= sc <= SC_MAX_A):
        reasons.append("sc_out_of_window")
    if oh >= OH_MAX_A:
        reasons.append("hydroxyl_oh_not_intact")
    if nonhydroxyl_min <= NONHYDROXYL_OC_MIN_A:
        reasons.append("nonhydroxyl_o_c_too_close")
    if mediator is None:
        reasons.append("no_mediating_water")
    if (d_iface_c is not None and d_iface_s is not None
            and max(d_iface_c, d_iface_s) > INTERFACE_MAX_A):
        reasons.append("not_near_interface")

    return {
        "step": record.step,
        "time_fs": record.time,
        "start_line": start_line,
        "sc_A": sc,
        "oh_A": oh,
        "o_nonhydroxyl_1_c_A": oc1,
        "o_nonhydroxyl_2_c_A": oc2,
        "o_nonhydroxyl_min_c_A": nonhydroxyl_min,
        "d5_A": None if mediator is None else mediator["d5_A"],
        "d3_A": None if mediator is None else mediator["d3_A"],
        "score_A": None if mediator is None else mediator["score_A"],
        "mediator": mediator,
        "frame_metrics_error": frame_metrics_error,
        "frame_metrics": metrics,
        "d_interface_c_A": d_iface_c,
        "d_interface_s_A": d_iface_s,
        "reasons": reasons,
        "accepted": not reasons,
    }


def find_mediator(atoms, box=CELL):
    """Lowest max(d5, d3) intact water bridging H393 -> Ow -> Hw -> O386.

    d5 = Ow...H393, d3 = O386...Hw.  Both water O-H bonds must be < 1.3 A, so
    a water that has already transferred a proton is not accepted as the intact
    three-atom mediator.
    """
    h_hydroxyl = atoms[IDX_H_HYDROXYL][1:4]
    o_carbonyl = atoms[IDX_O_CARBONYL][1:4]
    best = None
    for water in range(WATER_COUNT):
        oxygen_index = 3 * water
        h1_index = oxygen_index + 1
        h2_index = oxygen_index + 2
        oxygen = atoms[oxygen_index][1:4]
        if periodic_distance(oxygen, atoms[h1_index][1:4], box) >= MEDIATOR_OH_MAX_A:
            continue
        if periodic_distance(oxygen, atoms[h2_index][1:4], box) >= MEDIATOR_OH_MAX_A:
            continue
        d5 = periodic_distance(oxygen, h_hydroxyl, box)
        if d5 > BRIDGE_MAX_A:
            continue
        for h_index, label in ((h1_index, "H1"), (h2_index, "H2")):
            d3 = periodic_distance(atoms[h_index][1:4], o_carbonyl, box)
            if d3 > BRIDGE_MAX_A:
                continue
            candidate = {
                "water_index_0based": water,
                "oxygen_source_index": oxygen_index + 1,
                "hydrogen_source_indices": (h1_index + 1, h2_index + 1),
                "donor_h_label": label,
                "donor_h_source_index": h_index + 1,
                "donor_h_local_index": LOCAL_MEDIATOR_H1 if label == "H1" else LOCAL_MEDIATOR_H2,
                "d5_A": d5,
                "d3_A": d3,
                "score_A": max(d5, d3),
                "o_h1_A": periodic_distance(oxygen, atoms[h1_index][1:4], box),
                "o_h2_A": periodic_distance(oxygen, atoms[h2_index][1:4], box),
            }
            if (best is None
                    or candidate["score_A"] < best["score_A"]
                    or (candidate["score_A"] == best["score_A"]
                        and candidate["oxygen_source_index"] < best["oxygen_source_index"])
                    or (candidate["score_A"] == best["score_A"]
                        and candidate["oxygen_source_index"] == best["oxygen_source_index"]
                        and candidate["donor_h_source_index"] < best["donor_h_source_index"])):
                best = candidate
    return best


def cluster_indices(mediator):
    """0-based source indices in the fixed 13-atom order."""
    return list(SOLUTE_SOURCE_INDICES) + [
        mediator["oxygen_source_index"] - 1,
        mediator["hydrogen_source_indices"][0] - 1,
        mediator["hydrogen_source_indices"][1] - 1,
    ]


def cluster_atom_map(mediator):
    """Local -> source mapping, order preserving (no reordering)."""
    indices = cluster_indices(mediator)
    roles = list(SOLUTE_ROLES) + [
        "mediator_water_O", "mediator_water_H1", "mediator_water_H2"]
    elements = list(SOLUTE_ELEMENTS) + ["O", "H", "H"]
    mapping = []
    for local, (role, element, source_index) in enumerate(zip(roles, elements, indices)):
        mapping.append({
            "local": local,
            "local_1based": local + 1,
            "source_index_0based": source_index,
            "source_index_1based": source_index + 1,
            "role": role,
            "element": element,
            "mediator_donor_to_carbonyl": (
                local == mediator["donor_h_local_index"]),
        })
    return mapping


def extract_cluster(atoms, mediator):
    indices = cluster_indices(mediator)
    return [(atoms[index][0],) + tuple(float(value) for value in atoms[index][1:4])
            for index in indices]


def composition(atoms):
    return dict(Counter(atom[0] for atom in atoms))


def geometry_metrics(atoms, donor_local):
    """RC and water-bridge distances in the 13-atom cluster coordinate frame."""
    c = atoms[LOCAL_C][1:4]
    o_carbonyl = atoms[LOCAL_O_CARBONYL][1:4]
    s = atoms[LOCAL_S][1:4]
    o_hydroxyl = atoms[LOCAL_O_HYDROXYL][1:4]
    h_hydroxyl = atoms[LOCAL_H_HYDROXYL][1:4]
    o_non_1 = atoms[LOCAL_O_NONHYDROXYL_1][1:4]
    o_non_2 = atoms[LOCAL_O_NONHYDROXYL_2][1:4]
    mediator_o = atoms[LOCAL_MEDIATOR_O][1:4]
    donor_h = atoms[donor_local][1:4]
    return {
        "sc_A": periodic_distance(s, c),
        "oh_A": periodic_distance(o_hydroxyl, h_hydroxyl),
        "o_nonhydroxyl_1_c_A": periodic_distance(o_non_1, c),
        "o_nonhydroxyl_2_c_A": periodic_distance(o_non_2, c),
        "o_nonhydroxyl_min_c_A": min(periodic_distance(o_non_1, c),
                                     periodic_distance(o_non_2, c)),
        "d5_A": periodic_distance(mediator_o, h_hydroxyl),
        "d3_A": periodic_distance(o_carbonyl, donor_h),
        "bridge_max_A": max(periodic_distance(mediator_o, h_hydroxyl),
                            periodic_distance(o_carbonyl, donor_h)),
        "mediator_o_h1_A": periodic_distance(atoms[LOCAL_MEDIATOR_O][1:4],
                                             atoms[LOCAL_MEDIATOR_H1][1:4]),
        "mediator_o_h2_A": periodic_distance(atoms[LOCAL_MEDIATOR_O][1:4],
                                             atoms[LOCAL_MEDIATOR_H2][1:4]),
    }


def numpy_version():
    try:
        import numpy as np
    except ImportError:
        raise RcError("numpy is required for the Kabsch alignment")
    return str(np.__version__)


def kabsch_align(mobile, target):
    """Rigid Kabsch (SVD), returning (aligned, R, rmsd, det, mobile_centroid,
    target_centroid).  No scaling, no atom reordering."""
    try:
        import numpy as np
    except ImportError:
        raise RcError("numpy is required for the Kabsch alignment")
    points = np.asarray(mobile, dtype=float)
    reference = np.asarray(target, dtype=float)
    if points.shape != reference.shape or points.ndim != 2 or points.shape[1] != 3:
        raise RcError("Kabsch inputs must be matching N x 3 arrays")
    mobile_centroid = points.mean(axis=0)
    target_centroid = reference.mean(axis=0)
    centred = points - mobile_centroid
    reference_centred = reference - target_centroid
    covariance = centred.T.dot(reference_centred)
    u, _, vt = np.linalg.svd(covariance)
    determinant = float(np.sign(np.linalg.det(vt.T.dot(u.T))))
    if determinant == 0.0:
        raise RcError("degenerate Kabsch covariance (det = 0)")
    correction = np.diag([1.0, 1.0, determinant])
    rotation = vt.T.dot(correction).dot(u.T)
    aligned = centred.dot(rotation.T) + target_centroid
    difference = aligned - reference
    rmsd = float(np.sqrt(np.mean(np.sum(difference * difference, axis=1))))
    return aligned, rotation, rmsd, determinant, mobile_centroid, target_centroid


def apply_rigid(points, rotation, mobile_centroid, target_centroid):
    try:
        import numpy as np
    except ImportError:
        raise RcError("numpy is required for the Kabsch alignment")
    array = np.asarray(points, dtype=float)
    return (array - mobile_centroid).dot(rotation.T) + target_centroid


def min_distance_between_sets(atoms, first_indices, second_indices):
    best = None
    for first in first_indices:
        for second in second_indices:
            distance = periodic_distance(atoms[first][1:4], atoms[second][1:4])
            if best is None or distance < best["distance_A"]:
                best = {
                    "distance_A": distance,
                    "first_index_1based": first + 1,
                    "second_index_1based": second + 1,
                    "first_element": atoms[first][0],
                    "second_element": atoms[second][0],
                }
    return best


def search_clash_relief(positions, reference_atoms, cluster_indices,
                        max_shift_A=1.5, step_A=0.05,
                        min_contact_A=EMBED_MIN_EXTERNAL_A):
    """Bounded rigid translation that clears a cluster-external steric clash.

    The gas-phase-optimised 13-atom cluster is embedded by a Kabsch transform;
    after optimisation one of its atoms can sit inside a neighbouring solvent
    molecule (the 2026-10-28 stage-2 run landed 0.936 A from an external H).
    A rigid translation of the whole cluster along the clash axis is the minimal
    embedding refinement: the cluster's internal geometry is untouched and the
    >1.4 A contact gate, the unreacted-RC gates, the 2.5 A bridge gate and the
    5 A interface gate all stay exactly as pre-registered.  The search is
    bounded; if no shift <= max_shift_A clears the clash the caller keeps the
    unshifted placement and the normal gate rejects it.

    positions: list of [x, y, z] for the cluster, already Kabsch-aligned.
    reference_atoms: full frame, each (element, x, y, z).
    cluster_indices: 0-based indices of the cluster atoms in the frame.
    Returns (positions_maybe_shifted, info_dict).
    """
    cluster_set = set(cluster_indices)
    external = [i for i in range(len(reference_atoms)) if i not in cluster_set]

    def closest(points):
        best = None
        for local, pos in enumerate(points):
            for ext in external:
                ref = reference_atoms[ext]
                distance = periodic_distance(pos, ref[1:4])
                if best is None or distance < best["distance_A"]:
                    best = {"distance_A": distance, "local0": local,
                            "external0": ext}
        return best

    before = closest(positions)
    info = {
        "applied": False,
        "max_shift_A": max_shift_A,
        "step_A": step_A,
        "min_contact_gate_A": min_contact_A,
        "min_before_A": None if before is None else before["distance_A"],
    }
    if before is None or before["distance_A"] > min_contact_A:
        return positions, info

    # Direction candidates: the clash axis first, then the 26 unit directions of
    # a {-1,0,1}^3 stencil.  All are rigid translations of the whole cluster, so
    # the internal geometry and every gate/threshold stay identical.
    local = before["local0"]
    ref = reference_atoms[before["external0"]]
    delta = [0.0, 0.0, 0.0]
    for axis in range(3):
        d = positions[local][axis] - ref[axis + 1]
        box = CELL[axis]
        d -= box * round(d / box)
        delta[axis] = d
    norm = (delta[0] ** 2 + delta[1] ** 2 + delta[2] ** 2) ** 0.5
    directions = []
    if norm > 0.0:
        directions.append(("clash_axis", [value / norm for value in delta]))
    for dx in (-1.0, 0.0, 1.0):
        for dy in (-1.0, 0.0, 1.0):
            for dz in (-1.0, 0.0, 1.0):
                if dx == 0.0 and dy == 0.0 and dz == 0.0:
                    continue
                length = (dx ** 2 + dy ** 2 + dz ** 2) ** 0.5
                directions.append(("stencil", [dx / length, dy / length,
                                               dz / length]))

    shift = step_A
    while shift <= max_shift_A + 1e-9:
        for label, unit in directions:
            moved = [[pos[axis] + unit[axis] * shift for axis in range(3)]
                     for pos in positions]
            after = closest(moved)
            if after is not None and after["distance_A"] > min_contact_A:
                info.update({
                    "applied": True,
                    "translation_A": [value * shift for value in unit],
                    "magnitude_A": shift,
                    "min_after_A": after["distance_A"],
                    "clashing_local0": local,
                    "clashing_external0": before["external0"],
                    "direction_kind": label,
                    "direction_note": ("rigid translation of all cluster atoms; "
                                       "internal geometry, gates and thresholds "
                                       "unchanged"),
                })
                return moved, info
        shift += step_A
    return positions, info


def parse_orca_input(text):
    match = XYZFILE_RE.search(text)
    if match is None:
        raise RcError("ORCA input has no '* xyzfile charge multiplicity file' line")
    return {
        "charge": int(match.group(1)),
        "multiplicity": int(match.group(2)),
        "xyz_name": match.group(3),
    }


def parse_orca_output(text):
    """Parse termination, final energy, final geometry and all frequencies.

    The final geometry is taken from the block after the last
    'FINAL ENERGY EVALUATION AT THE STATIONARY POINT' marker; falling back to
    an arbitrary coordinate block is deliberately not supported.
    """
    energies = ENERGY_RE.findall(text)
    if not energies:
        raise RcError("ORCA output has no FINAL SINGLE POINT ENERGY")
    final_energy = float(energies[-1])

    anchor = text.rfind("FINAL ENERGY EVALUATION AT THE STATIONARY POINT")
    if anchor < 0:
        raise RcError(
            "ORCA output has no FINAL ENERGY EVALUATION AT THE STATIONARY POINT")
    coord_match = COORD_BLOCK_RE.search(text, anchor)
    if coord_match is None:
        raise RcError("no Cartesian coordinate block after the stationary-point marker")
    atoms = []
    for line in coord_match.group(1).splitlines():
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) != 4:
            raise RcError("bad final-geometry row in ORCA output: {!r}".format(line))
        if not ELEMENT_RE.match(fields[0]):
            raise RcError("bad element in ORCA final geometry: {!r}".format(fields[0]))
        try:
            atoms.append((fields[0], float(fields[1]), float(fields[2]),
                          float(fields[3])))
        except ValueError:
            raise RcError("non-numeric final-geometry row: {!r}".format(line))
    if not atoms:
        raise RcError("empty final-geometry block")

    freq_start = text.rfind("VIBRATIONAL FREQUENCIES")
    if freq_start < 0:
        raise RcError("ORCA output has no VIBRATIONAL FREQUENCIES block")
    freq_end = text.find("NORMAL MODES", freq_start)
    freq_section = text[freq_start:freq_end if freq_end > 0 else len(text)]
    pairs = [(int(match.group(1)), float(match.group(2)))
             for match in FREQ_LINE_RE.finditer(freq_section)]
    if not pairs:
        raise RcError("VIBRATIONAL FREQUENCIES block has no modes")
    mode_ids = [mode for mode, _ in pairs]
    if mode_ids != list(range(len(pairs))):
        raise RcError("frequency mode indices are not contiguous from 0")
    frequencies = [value for _, value in pairs]

    imaginary_report = re.findall(
        r"Total number of imaginary perturbations\s*\.\.\.\s*(\d+)", text)
    return {
        "terminated_normally": "ORCA TERMINATED NORMALLY" in text,
        "hurray_converged": ("HURRAY" in text
                             and "THE OPTIMIZATION HAS CONVERGED" in text),
        "final_single_point_energy_Eh": final_energy,
        "final_energy_count": len(energies),
        "final_geometry": atoms,
        "frequencies_cm1": frequencies,
        "reported_imaginary_perturbations": (
            int(imaginary_report[-1]) if imaginary_report else None),
        "has_stationary_point_marker": True,
    }


def frequency_gate(frequencies, natoms, threshold=IMAGINARY_FAIL_BELOW_CM1,
                   reported_imaginary=None):
    """No frequency below -10 cm-1 passes.  [-10, 0) is recorded noise.

    ``reported_imaginary`` is ORCA's own
    ``Total number of imaginary perturbations`` count.  A positive count fails
    the gate even when no negative frequency line was parsed, so a parse change
    cannot turn a real imaginary mode into a silent pass.
    """
    expected_modes = 3 * natoms
    negative = [value for value in frequencies if value < -1.0e-12]
    fail = [value for value in frequencies if value < threshold]
    noise = [value for value in frequencies if -1.0e-12 > value >= threshold]
    reported_fail = reported_imaginary is not None and reported_imaginary > 0
    return {
        "expected_modes": expected_modes,
        "n_modes": len(frequencies),
        "mode_count_ok": len(frequencies) == expected_modes,
        "threshold_cm1": threshold,
        "frequencies_cm1": frequencies,
        "negative_frequencies_cm1": negative,
        "numerical_noise_frequencies_cm1": noise,
        "imaginary_below_threshold_cm1": fail,
        "reported_imaginary_perturbations": reported_imaginary,
        "reported_imaginary_fails": reported_fail,
        "passes": (len(frequencies) == expected_modes and not fail
                   and not reported_fail),
        "policy": ("any frequency < -10 cm-1 fails; frequencies in [-10, 0) are "
                   "reported as numerical noise and are never silently ignored; "
                   "a positive ORCA 'Total number of imaginary perturbations' "
                   "also fails"),
    }

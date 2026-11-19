#!/usr/bin/env python3
"""Export Li 2025 pathway-3 (SO3^2- + HCHO) R snapshots and build CP2K/PLUMED inputs.

Provenance
----------
J. Li et al., Nature Communications 16, 5187 (2025), DOI 10.1038/s41467-025-59712-3.
Zenodo archive DOI 10.5281/zenodo.15253851 (CC BY 4.0). The workbook lives at
results/li2025_zenodo/source-data.xlsx and its MD5 is pinned below.

Only sheet "Figure 6b" is read, and only the two reactant blocks:
  R_bulk  -> compound columns A:D  (element A, xyz B:D)
  R_inter -> compound columns P:S  (element P, xyz Q:S)
Rows 2..297 hold 296 atoms. Coordinates are copied verbatim: no wrapping,
translation, rotation, atom deletion or reordering. The stored values are
six-decimal literals, so "%.6f" reproduces them exactly (checked below).

The public snapshot is 96 H2O + SO3^2- + HCHO with no K+ / Na+. This is therefore
a periodic q = -2 model. It is NOT the project's +2 K+ neutral 128-water
PBE-D3 model and must never be merged with it, nor quoted as the same system.

This script also refuses to emit anything when the workbook, the composition,
the atom order, the solute connectivity or the contact geometry fails a gate.
Exact duplicate atoms and unparseable rows are hard errors: we do not "fix"
coordinates to chase a number.

Outputs (all inside campaign/li2025_direct/):
  inputs/bulk_r1.xyz,      inputs/bulk_r1.inp,      inputs/bulk_r1.plumed
  inputs/interface_r1.xyz, inputs/interface_r1.inp, inputs/interface_r1.plumed
  manifest.json

Run: python3 prepare_li2025_direct.py
"""
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"
SOURCE = ROOT / "results" / "li2025_zenodo" / "source-data.xlsx"

SOURCE_MD5 = "53a0f1915b8f15fe1cd5d56742d16820"
PAPER_DOI = "10.1038/s41467-025-59712-3"
ZENODO_DOI = "10.5281/zenodo.15253851"
SHEET = "Figure 6b"
ROW_FIRST = 2
ROW_LAST = 297
ATOM_COUNT = 296

BLOCKS = [("R_bulk", 1), ("R_inter", 16)]
BLOCK_COLUMN_LETTERS = {
    "R_bulk": {"element": "A", "xyz": "B:D", "cell": "A2:D297"},
    "R_inter": {"element": "P", "xyz": "Q:S", "cell": "P2:S297"},
}

EXPECTED_COMPOSITION = {"H": 194, "O": 100, "C": 1, "S": 1}
EXPECTED_ORDER = ["H", "H", "O"] * 96 + ["C", "O", "H", "H", "S", "O", "O", "O"]
SOLUTE = {"C": 289, "O_carbonyl": 290, "H_formyl_1": 291, "H_formyl_2": 292,
          "S": 293, "O_sulfonate_1": 294, "O_sulfonate_2": 295, "O_sulfonate_3": 296}

# 1 Ha = 627.509 kcal/mol is the factor behind both requested HEIGHT literals.
HA_TO_KCAL = 627.509
WALL_AT = 4.0
WALL_KAPPA = 200.0

CASES = [
    {
        "name": "bulk_r1",
        "point": "R_bulk",
        "environment": "bulk",
        "column": 1,
        "box": (14.2296, 14.2296, 14.2296),
        "height_ha": 1.0e-4,
        "height_kcal": 0.0627509,
    },
    {
        "name": "interface_r1",
        "point": "R_inter",
        "environment": "interface",
        "column": 16,
        "box": (15.0, 15.0, 30.0),
        "height_ha": 5.0e-5,
        "height_kcal": 0.03137545,
    },
]

CP2K_INPUT = """&GLOBAL
  PROJECT {project}
  RUN_TYPE MD
  PRINT_LEVEL MEDIUM
&END GLOBAL
&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME /public/home/temp/cp2k-9.1/data/BASIS_MOLOPT
    POTENTIAL_FILE_NAME /public/home/temp/cp2k-9.1/data/GTH_POTENTIALS
    CHARGE -2
    MULTIPLICITY 1
    &MGRID
      CUTOFF 300
      REL_CUTOFF 40
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-12
    &END QS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-6
      MAX_SCF 50
      &OT
        MINIMIZER DIIS
      &END OT
      &OUTER_SCF
        MAX_SCF 20
        EPS_SCF 1.0E-6
      &END OUTER_SCF
    &END SCF
    &XC
      &XC_FUNCTIONAL BLYP
      &END XC_FUNCTIONAL
      &VDW_POTENTIAL
        POTENTIAL_TYPE PAIR_POTENTIAL
        &PAIR_POTENTIAL
          TYPE DFTD3(BJ)
          PARAMETER_FILE_NAME /public/home/temp/cp2k-9.1/data/dftd3.dat
          REFERENCE_FUNCTIONAL BLYP
        &END PAIR_POTENTIAL
      &END VDW_POTENTIAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC {a} {b} {c}
      PERIODIC XYZ
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME {project}.xyz
      COORD_FILE_FORMAT XYZ
    &END TOPOLOGY
    &KIND H
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-BLYP-q1
    &END KIND
    &KIND O
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-BLYP-q6
    &END KIND
    &KIND C
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-BLYP-q4
    &END KIND
    &KIND S
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-BLYP-q6
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
&MOTION
  &FREE_ENERGY
    &METADYN
      USE_PLUMED .TRUE.
      PLUMED_INPUT_FILE ./{project}.plumed
    &END METADYN
  &END FREE_ENERGY
  &MD
    ENSEMBLE NVT
    STEPS 10000
    TIMESTEP 1.0
    TEMPERATURE 298.15
    &THERMOSTAT
      TYPE CSVR
      &CSVR
        TIMECON 100
      &END CSVR
    &END THERMOSTAT
  &END MD
  &PRINT
    &TRAJECTORY
      FORMAT XYZ
      &EACH
        MD 10
      &END EACH
    &END TRAJECTORY
    &RESTART
      &EACH
        MD 100
      &END EACH
    &END RESTART
  &END PRINT
&END MOTION
"""

PLUMED_INPUT = """# Li 2025 Fig. 6b {point}: SO3^2- + HCHO direct C-S addition, single d(S-C) CV.
# First production replica ({project}); ordinary (non-well-tempered) metaD.
# PACE and PRINT STRIDE are MD steps (1 fs each), so 30 steps = 30 fs, 10 steps = 10 fs.
UNITS LENGTH=A ENERGY=kcal/mol TIME=fs
d_sc: DISTANCE ATOMS={s},{c}
metad: METAD ARG=d_sc SIGMA=0.1 HEIGHT={height} PACE=30 FILE=HILLS_{project}
PRINT ARG=d_sc,metad.bias STRIDE=10 FILE=COLVAR_{project}
uw_sc: UPPER_WALLS ARG=d_sc AT={wall_at} KAPPA={wall_kappa}
"""


def md5(path):
    digest = hashlib.md5()
    with open(str(path), "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def distance(a, b, box=None):
    d = [a[i] - b[i] for i in range(3)]
    if box is not None:
        d = [d[i] - box[i] * round(d[i] / box[i]) for i in range(3)]
    return math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])


def wrapped(point, box):
    return [p - box[i] * math.floor(p / box[i]) for i, p in enumerate(point)]


def read_block(worksheet, column):
    """Read one Figure 6b block. column is the 1-based element column."""
    rows = list(worksheet.iter_rows(min_row=ROW_FIRST, max_row=ROW_LAST, values_only=True))
    if len(rows) != ATOM_COUNT:
        raise ValueError("expected %d rows, got %d" % (ATOM_COUNT, len(rows)))
    elements = []
    coords = []
    for offset, row in enumerate(rows):
        row_number = ROW_FIRST + offset
        element = row[column - 1]
        if not isinstance(element, str) or element.strip() not in ("H", "C", "O", "S"):
            raise ValueError("row %d: unparseable element %r" % (row_number, element))
        xyz = row[column:column + 3]
        values = []
        for value in xyz:
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise ValueError("row %d: unparseable coordinate %r" % (row_number, value))
            if not math.isfinite(float(value)):
                raise ValueError("row %d: non-finite coordinate %r" % (row_number, value))
            values.append(float(value))
        elements.append(element.strip())
        coords.append(values)
    return elements, coords


def molecules(elements):
    """Fixed grouping for this known order: 96 waters, then HCHO, then SO3."""
    waters = []
    for start in range(0, 288, 3):
        waters.append([start, start + 1, start + 2])
    return {
        "water": waters,
        "hcho": [288, 289, 290, 291],
        "so3": [292, 293, 294, 295],
    }


def verify_structure(elements, coords, box, label):
    """Hard gates. Raises ValueError instead of editing the structure."""
    metrics = {}
    if len(elements) != ATOM_COUNT:
        raise ValueError("%s: expected %d atoms, got %d" % (label, ATOM_COUNT, len(elements)))
    composition = dict(Counter(elements))
    if composition != EXPECTED_COMPOSITION:
        raise ValueError("%s: composition %s != %s" % (label, composition, EXPECTED_COMPOSITION))
    if elements != EXPECTED_ORDER:
        raise ValueError("%s: atom order differs from the expected "
                         "96x(H,H,O) + C,O,H,H,S,O,O,O order" % label)

    # Exact duplicate atoms, including duplicates modulo the simulation cell.
    closest = (1.0e9, None)
    for i in range(ATOM_COUNT):
        for j in range(i + 1, ATOM_COUNT):
            d = distance(coords[i], coords[j], box)
            if d < closest[0]:
                closest = (d, (i + 1, j + 1, elements[i], elements[j]))
    if closest[0] < 1.0e-6:
        raise ValueError("%s: exact duplicate atoms at 1-based indices %s (%.3e A)"
                         % (label, closest[1], closest[0]))
    metrics["min_atom_pair_distance_angstrom"] = closest[0]
    metrics["min_atom_pair_1based"] = list(closest[1])

    groups = molecules(elements)
    for water in groups["water"]:
        oxygen = water[2]
        hs = [h for h in water if elements[h] == "H"]
        if len(hs) != 2 or any(distance(coords[oxygen], coords[h], box) > 1.30 for h in hs):
            raise ValueError("%s: water at 1-based O%d is not an intact H2O" % (label, oxygen + 1))
    s_index = 292
    s_oxygen = [o for o in groups["so3"] if elements[o] == "O"]
    if len(s_oxygen) != 3 or any(distance(coords[s_index], coords[o], box) > 1.75 for o in s_oxygen):
        raise ValueError("%s: SO3 is not an intact S + 3 O unit" % label)
    c_index = 288
    c_oxygen = [o for o in groups["hcho"] if elements[o] == "O"]
    c_hydrogen = [h for h in groups["hcho"] if elements[h] == "H"]
    if (len(c_oxygen) != 1 or distance(coords[c_index], coords[c_oxygen[0]], box) > 1.45
            or len(c_hydrogen) != 2
            or any(distance(coords[c_index], coords[h], box) > 1.25 for h in c_hydrogen)):
        raise ValueError("%s: HCHO is not an intact C=O + 2 H unit" % label)

    # Shortest intermolecular contact (intramolecular bonds are excluded by group).
    mol_list = [list(water) for water in groups["water"]]
    mol_labels = ["H2O"] * len(mol_list)
    mol_list.append(groups["hcho"])
    mol_labels.append("HCHO")
    mol_list.append(groups["so3"])
    mol_labels.append("SO3")
    best = (1.0e9, None)
    for a in range(len(mol_list)):
        for b in range(a + 1, len(mol_list)):
            for i in mol_list[a]:
                for j in mol_list[b]:
                    d = distance(coords[i], coords[j], box)
                    if d < best[0]:
                        best = (d, (i + 1, mol_labels[a], j + 1, mol_labels[b],
                                    elements[i], elements[j]))
    if best[0] < 1.20:
        raise ValueError("%s: severe overlap, intermolecular contact %.3f A between %s"
                         % (label, best[0], best[1]))
    metrics["min_intermolecular_contact_angstrom"] = best[0]
    metrics["min_intermolecular_contact_1based"] = list(best[1])

    metrics["sc_distance_angstrom"] = distance(coords[292], coords[288], box)
    metrics["composition"] = composition
    return metrics


def film_geometry(coords, elements, box):
    """Largest O-z gap of the interface snapshot: film thickness and vacuum gap."""
    oxygen_z = sorted(wrapped(coords[i], box)[2] for i in range(288) if elements[i] == "O")
    gaps = [oxygen_z[i + 1] - oxygen_z[i] for i in range(len(oxygen_z) - 1)]
    index = max(range(len(gaps)), key=lambda k: gaps[k])
    bottom, top = oxygen_z[index], oxygen_z[index + 1]
    thickness = box[2] - gaps[index]
    return {
        "water_oxygen_z_min_angstrom": oxygen_z[0],
        "water_oxygen_z_max_angstrom": oxygen_z[-1],
        "largest_water_oxygen_gap_angstrom": gaps[index],
        "film_thickness_angstrom": thickness,
        "film_top_surface_z_angstrom": bottom,
        "film_bottom_surface_z_angstrom": top,
        "vacuum_gap_angstrom": gaps[index],
        "solute_wrapped_z_angstrom": {
            "S": wrapped(coords[292], box)[2],
            "C": wrapped(coords[288], box)[2],
        },
        "carbon_depth_below_top_surface_angstrom": bottom - wrapped(coords[288], box)[2],
        "water_density_g_cm3_in_film": (
            96 * 18.015 / 6.02214076e23
            / (box[0] * 1.0e-8 * box[1] * 1.0e-8 * thickness * 1.0e-8) / 1.0
        ),
    }


def xyz_text(elements, coords, comment):
    lines = [str(len(elements)), comment]
    for element, point in zip(elements, coords):
        lines.append("%-2s %12.6f %12.6f %12.6f" % (element, point[0], point[1], point[2]))
    return "\n".join(lines) + "\n"


def build_case(elements, coords, case, metrics):
    name = case["name"]
    box = case["box"]
    project = name
    xyz_path = INPUTS / (name + ".xyz")
    inp_path = INPUTS / (name + ".inp")
    plumed_path = INPUTS / (name + ".plumed")

    comment = ("Li2025 Fig6b %s | 296 atoms = 96 H2O + SO3^2- + HCHO | q=-2, no K/Na"
               " | coords verbatim from Zenodo source-data.xlsx" % case["point"])
    xyz_path.write_text(xyz_text(elements, coords, comment), encoding="utf-8")
    inp_path.write_text(CP2K_INPUT.format(project=project, a=box[0], b=box[1], c=box[2]),
                        encoding="utf-8")
    plumed_path.write_text(PLUMED_INPUT.format(
        point=case["point"], project=project, s=SOLUTE["S"], c=SOLUTE["C"],
        height=case["height_kcal"], wall_at=WALL_AT, wall_kappa=WALL_KAPPA), encoding="utf-8")

    record = {
        "name": name,
        "point": case["point"],
        "environment": case["environment"],
        "xyz": "inputs/%s.xyz" % name,
        "inp": "inputs/%s.inp" % name,
        "plumed": "inputs/%s.plumed" % name,
        "xyz_sha256": sha256(xyz_path),
        "inp_sha256": sha256(inp_path),
        "plumed_sha256": sha256(plumed_path),
        "atom_count": ATOM_COUNT,
        "composition": metrics["composition"],
        "box_angstrom": list(box),
        "sc_distance_angstrom": metrics["sc_distance_angstrom"],
        "min_intermolecular_contact_angstrom": metrics["min_intermolecular_contact_angstrom"],
        "min_intermolecular_contact_atoms_1based": metrics["min_intermolecular_contact_1based"],
        "min_intermolecular_contact_metric": ("minimum-image distance, intramolecular pairs "
                                              "excluded; hard gate at 1.20 A"),
        "min_atom_pair_distance_including_bonds_angstrom":
            metrics["min_atom_pair_distance_angstrom"],
        "min_atom_pair_including_bonds_1based": metrics["min_atom_pair_1based"],
        "max_coordinate_excursion_angstrom": max(
            [0.0] + [max(0.0, -point[i], point[i] - box[i])
                     for point in coords for i in range(3)]),
        "solute_positions_angstrom_raw": {
            "S": coords[292], "C": coords[288],
            "O_carbonyl": coords[289],
        },
        "metad": {
            "cv": "d(S-C) single distance",
            "sigma_angstrom": 0.1,
            "pace_md_steps": 30,
            "pace_fs_at_1fs_timestep": 30,
            "print_stride_md_steps": 10,
            "print_stride_fs_at_1fs_timestep": 10,
            "height_ha": case["height_ha"],
            "height_kcal_mol": case["height_kcal"],
            "biasfactor": None,
            "note": ("ordinary metaD; no BIASFACTOR and no GRID/WT block. PLUMED PACE and "
                     "PRINT STRIDE count MD steps, which equal fs here because TIMESTEP=1.0."),
        },
        "wall": {
            "type": "UPPER_WALLS",
            "at_angstrom": WALL_AT,
            "kappa_kcal_mol_angstrom2": WALL_KAPPA,
            "initial_margin_angstrom": WALL_AT - metrics["sc_distance_angstrom"],
            "initial_wall_energy_kcal_mol": 0.0,
            "status": ("AT = 4.0 A matches SI Table 3; KAPPA = 200 is a project assumption "
                       "because the wall strength is not reported"),
        },
        "cp2k": {
            "charge": -2,
            "multiplicity": 1,
            "ensemble": "NVT",
            "thermostat": "CSVR",
            "temperature_k": 298.15,
            "timecon_fs": 100,
            "timestep_fs": 1.0,
            "steps_initial_segment": 10000,
            "trajectory_every_fs": 10,
            "restart_every_fs": 100,
        },
    }
    if case["environment"] == "interface":
        record["interface_geometry"] = film_geometry(coords, elements, box)
    else:
        volume_cm3 = (box[0] * box[1] * box[2]) * 1.0e-24
        record["bulk_geometry"] = {
            "water_density_g_cm3": 96 * 18.015 / 6.02214076e23 / volume_cm3,
            "water_oxygen_z_min_angstrom": min(
                wrapped(coords[i], box)[2] for i in range(288) if elements[i] == "O"),
            "water_oxygen_z_max_angstrom": max(
                wrapped(coords[i], box)[2] for i in range(288) if elements[i] == "O"),
        }
    return record


def main():
    INPUTS.mkdir(parents=True, exist_ok=True)
    actual_md5 = md5(SOURCE)
    if actual_md5 != SOURCE_MD5:
        raise SystemExit("workbook MD5 mismatch: %s != %s" % (actual_md5, SOURCE_MD5))

    from openpyxl import load_workbook

    workbook = load_workbook(str(SOURCE), read_only=True, data_only=True)
    worksheet = workbook[SHEET]
    header = [cell.value for cell in next(worksheet.iter_rows(min_row=1, max_row=1))]
    for label, column in BLOCKS:
        if header[column - 1] != label:
            raise SystemExit("sheet %s column %d header %r != %r"
                             % (SHEET, column, header[column - 1], label))

    payload = {}
    for label, column in BLOCKS:
        elements, coords = read_block(worksheet, column)
        box = next(case["box"] for case in CASES if case["point"] == label)
        metrics = verify_structure(elements, coords, box, label)
        payload[label] = (elements, coords, metrics)
    workbook.close()

    records = []
    for case in CASES:
        elements, coords, metrics = payload[case["point"]]
        records.append(build_case(elements, coords, case, metrics))

    manifest = {
        "session": "li2025-direct-pathway3",
        "purpose": ("First production replicas of the Li 2025 pathway-3 (SO3^2- + HCHO) "
                    "direct C-S addition, started from the published bulk and interface "
                    "reactant snapshots."),
        "scope": {
            "what_this_is": [
                "Two independent production workdirs (bulk_r1_fix1, interface_r1_fix1, 64 MPI ranks each).",
                "R snapshots only; the TS and P blocks of Figure 6b are not exported.",
                "Ordinary single-CV metaD on d(S-C), 10 ps initial segment per replica.",
            ],
            "what_this_is_not": [
                "NOT Li's three-replica free-energy statistics; this is replica 1 of 2 workdirs.",
                "NOT a Liu 2026 five-distance-CV calculation and must not be compared with it.",
                "No FES/PMF exists yet. HILLS/COLVAR appear only after a production job runs, "
                "and 10 ps at these heights is not a converged free-energy surface.",
                "NOT the project's +2 K+ neutral model: this public snapshot has no counterions.",
            ],
            "structure_conclusion_status": ("No structural conclusion is drawn here. Any later "
                                            "C-S bonding claim must pass the three-layer gate "
                                            "in AGENTS.md section 7."),
        },
        "source": {
            "paper": "J. Li et al., Nature Communications 16, 5187 (2025)",
            "doi": PAPER_DOI,
            "zenodo_doi": ZENODO_DOI,
            "license": "CC BY 4.0",
            "workbook": "results/li2025_zenodo/source-data.xlsx",
            "workbook_md5": actual_md5,
            "workbook_sha256": sha256(SOURCE),
            "sheet": SHEET,
            "rows": "%d:%d" % (ROW_FIRST, ROW_LAST),
            "atom_count_per_block": ATOM_COUNT,
            "blocks_exported": BLOCK_COLUMN_LETTERS,
            "extraction": ("Element and xyz read per row from the block columns; nothing is "
                           "wrapped, translated, rotated or reordered. Stored values are "
                           "six-decimal literals (max |v-round(v,6)| = 1.8e-15), so %.6f "
                           "is an exact reproduction."),
        },
        "system": {
            "stoichiometry": "96 H2O + SO3^2- + HCHO",
            "composition": EXPECTED_COMPOSITION,
            "net_charge": -2,
            "multiplicity": 1,
            "counterions": "none (no K+ or Na+ in this public snapshot)",
            "charge_model_note": (
                "q = -2 periodic model. The two sulfonate charges are not compensated by "
                "counterions, so CP2K applies a homogeneous neutralizing background. This "
                "system is not the project's +2 K+ neutral model, and energies, densities or "
                "barriers from the two must not be merged or compared as if they were the "
                "same calculation."),
            "atom_order": ("1-288 = 96 waters in H,H,O order (O at 3,6,...,288); "
                           "289 C, 290 O, 291 H, 292 H (HCHO); 293 S, 294-296 O (SO3)"),
            "solute_indices_1based": SOLUTE,
        },
        "production": {
            "engine": "CP2K 9.1 (cp2k.psmp with PLUMED)",
            "binary": "/public/home/temp/cp2k-9.1/exe/Linux-x86-64-gfortran-mkl-plumed/cp2k.psmp",
            "plumed": "2.7.3 (shared build, /public/home/temp/apps/plumed-2.7.3-shared/lib)",
            "functional": "BLYP",
            "dispersion": ("DFTD3(BJ) with REFERENCE_FUNCTIONAL BLYP (Becke-Johnson damping). "
                           "Plain 'TYPE DFTD3' would silently mean zero damping; see "
                           "docs/methods/mlp_plan.md lines 90-91."),
            "basis": "DZVP-MOLOPT-SR-GTH",
            "potentials": {"H": "GTH-BLYP-q1", "O": "GTH-BLYP-q6",
                           "C": "GTH-BLYP-q4", "S": "GTH-BLYP-q6"},
            "cutoff_ry": 300,
            "rel_cutoff_ry": 40,
            "scf": {"eps_scf": 1.0e-6, "ot_minimizer": "DIIS", "outer_scf": True},
            "ensemble": "NVT",
            "thermostat": "CSVR",
            "temperature_k": 298.15,
            "timecon_fs": 100,
            "timestep_fs": 1.0,
            "steps_initial_segment": 10000,
            "trajectory_every_fs": 10,
            "restart_every_fs": 100,
            "mpi_ranks": 64,
            "memory": "32G",
            "max_walltime_hours": 48,
            "partition": "compute",
            "mail_type": "END,FAIL",
            "assumptions": [
                "CSVR TIMECON = 100 fs is a project assumption; Li 2025 does not state it.",
                "UPPER_WALLS AT = 4.0 A matches SI Table 3; KAPPA = 200 is a project assumption.",
                "HEIGHT 1e-4 Ha (bulk) and 0.5e-4 Ha (interface), SIGMA = 0.1 A, PACE = 30 fs "
                "match the two pathway-3 rows in SI Table 3; both rows are printed as 3a, so "
                "the interface assignment of the 0.5e-4 Ha row follows the footnote scheme.",
                "Our implementation uses ordinary fixed-height hills; the SI does not establish "
                "whether the published calculations used ordinary or well-tempered metaD.",
                "1 Ha = 627.509 kcal/mol reproduces both requested HEIGHT literals.",
                "Only the first 10 ps segment is defined here; continuation is a separate step.",
            ],
        },
        "fes": {
            "status": "not_computed",
            "note": ("Nothing in this directory is a free-energy surface. HILLS/COLVAR are "
                     "raw metaD output and would need convergence and replica analysis before "
                     "any barrier is quoted."),
        },
        "cases": records,
    }
    (HERE / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    summary = {
        "workbook_md5": actual_md5,
        "cases": [
            {
                "name": record["name"],
                "sc_distance_angstrom": round(record["sc_distance_angstrom"], 4),
                "min_intermolecular_contact_angstrom":
                    round(record["min_intermolecular_contact_angstrom"], 4),
                "box_angstrom": record["box_angstrom"],
                "xyz_sha256": record["xyz_sha256"],
            }
            for record in records
        ],
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

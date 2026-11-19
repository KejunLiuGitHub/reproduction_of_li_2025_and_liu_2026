#!/usr/bin/env python3
"""D18 C-S coupling analysis: primary criterion = C-S bond only.

The user (2026-10-28 16:2x) changed the primary visible event to **C-S bond
formation**: "S-C formation is enough, proton transfer does not matter".  This
analyzer therefore answers a DIFFERENT question from ``analyze_unbiased.py``:

    criterion "cs_bond_only":  S(388)-C(384) < 2.0 A, sustained
                               (>= MIN_PERSISTENCE_FRAMES = 10 frames) to the
                               last frame; proton transfer is NOT required.

It is **not** the paper-comparable criterion and the JSON says so explicitly
(``criterion_is_paper_comparable: false``).  The original
``analyze_unbiased`` product definition (S-C < 2.0 A **and** O-H > 1.5 A) is
still computed and reported separately under ``paper_criterion`` so both numbers
can be quoted side by side.  Nothing about the old criterion is relaxed or
re-scored: it is simply no longer the primary event here.

What is imported (never re-typed) from ``analyze_unbiased``:
* the strict parser contract (comment ``i = <step>``/``time = <value>``),
  periodic-distance helper, frame dedup policy and the C/S/O-H element checks;
* the product thresholds ``S_C_FORMED_MAX_A`` / ``O_H_BROKEN_MIN_A``;
* the sustained-persistence rule ``MIN_PERSISTENCE_FRAMES`` and the
  3 ps sub-bar ``MIN_COVERAGE_FS``;
* the coverage gates (first step 0 / time 0, >=3000 fs span, >=300 unique
  frames, <=20 fs adjacent gap) and the temperature gates via
  ``coverage_block`` / ``temperature_block``.

In-progress trajectories
------------------------
The loader tolerates a trajectory that is still being written by CP2K: trailing
incomplete frames are skipped (recorded in ``trajectory_read``), and if the last
candidate frame cannot be parsed strictly it is treated as "still being written"
and dropped once (a second failure is a structural error -> exit 2).

Secondary observables (reported, never used as the primary gate)
----------------------------------------------------------------
O-H min/max/final, whether proton transfer happened (O-H > 1.5 A sustained),
the longest O-H > 1.5 A run, and the HMS diagnostic distances d2 (water O-H),
d3 (water donor H ... O386) and d5 (water O ... H393) as summaries, when the
atom map is available (seed manifest / strong-CS manifest provides it).

Verdicts / exit codes
---------------------
0  CS_COUPLING_SUSTAINED_BOND
4  CS_COUPLING_INITIAL_STATE_NOT_UNREACTED
4  CS_COUPLING_NO_SUSTAINED_BOND
4  CS_COUPLING_SUSTAINED_TEMPERATURE_FAIL
5  CS_COUPLING_COVERAGE_SHORT
2  structural / provenance error (bad trajectory, unknown biased flag, ...)

Python 3.6 compatible: the cluster interpreter is old.
"""
import argparse
import json
import tempfile
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import analyze_iface  # noqa: E402
import analyze_unbiased as ub  # noqa: E402  (strict parser + reused thresholds)
import iface_rc_common as common  # noqa: E402  (local cluster-role constants)

PROTOCOL = "D18_cs_coupling_analysis"
CRITERION = "cs_bond_only"
DEFAULT_JSON = "CS_COUPLING.json"
DEFAULT_TSV = "CS_COUPLING.tsv"

MANIFEST_CANDIDATES = (
    "IFACE_STRONG_CS_INPUT.json",     # D18 strong S-C steer
    "IFACE_STEERED_INPUT.json",       # D15 steered BOMD
    "IFACE_UNBIASED_INPUT.json",      # unbiased 3 ps BOMD
)
SEED_MANIFEST_CANDIDATES = ("IFACE_RC_SEED.json",)
TRAJECTORY_CANDIDATES = (
    "iface_strong_cs-pos-1.xyz",
    "iface_steered-pos-1.xyz",
    "iface_unbiased-pos-1.xyz",
)
ENERGY_CANDIDATES = (
    "iface_strong_cs-1.ener",
    "iface_steered-1.ener",
    "iface_unbiased-1.ener",
)
PROTOCOL_BIASED = {
    "D18_strong_cs_steer": True,
    "D15_steered_bomd": True,
}
PROTOCOL_UNBIASED_MARKERS = ("iface_unbiased",)

VERDICTS = {
    0: "CS_COUPLING_SUSTAINED_BOND",
    4: "CS_COUPLING_NO_SUSTAINED_BOND",
    5: "CS_COUPLING_COVERAGE_SHORT",
}
FORBIDDEN_OUTPUT_TOKENS = ("CRITERION_SATISFIED", "CRITERION_FAILED")

ROLE_KEYS = {
    "C": "HCHO_C",
    "O_carbonyl": "HCHO_O_carbonyl",
    "S": "HOSO2_S",
    "O_hydroxyl": "HOSO2_O_hydroxyl",
    "H_hydroxyl": "HOSO2_H_hydroxyl",
    "water_O": "mediator_water_O",
    "water_H1": "mediator_water_H1",
    "water_H2": "mediator_water_H2",
}


class StructureError(RuntimeError):
    """Structural/provenance error that maps to exit code 2."""


def sha256(path):
    import hashlib
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1048576), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path, label):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise StructureError("{} is missing or empty: {}".format(label, path))
    try:
        data = json.loads(path.read_text())
    except ValueError as exc:
        raise StructureError("{} is not valid JSON ({}): {}".format(
            label, exc, path))
    if not isinstance(data, dict):
        raise StructureError("{} is not a JSON object: {}".format(label, path))
    return data


def frame_prefix(lines):
    """Return (frame_start_end_pairs, partial_frame_info).

    A frame is header line (atom count) + comment line + N atom lines.  Any
    trailing block with fewer atom lines than its header declares is reported as
    a partial frame (the file is still being written) instead of a hard error.
    """
    frames = []
    partial = None
    index = 0
    total = len(lines)
    while index < total:
        stripped = lines[index].strip()
        if not stripped:
            index += 1
            continue
        try:
            natoms = int(stripped)
        except ValueError:
            raise StructureError(
                "line {}: XYZ frame header is not an atom count: {!r}".format(
                    index + 1, stripped))
        if natoms <= 0:
            raise StructureError(
                "line {}: atom count must be positive, got {}".format(
                    index + 1, natoms))
        if index + 2 > total:
            partial = {
                "start_line": index + 1,
                "expected_atoms": natoms,
                "atom_lines_present": 0,
                "reason": "file ends before the XYZ comment line",
            }
            break
        needed = index + 2 + natoms
        if needed > total:
            partial = {
                "start_line": index + 1,
                "expected_atoms": natoms,
                "atom_lines_present": max(0, total - (index + 2)),
                "reason": "file ends mid-frame",
            }
            break
        frames.append((index, needed))
        index = needed
    return frames, partial


def load_trajectory_in_progress(path, box):
    """Tolerant loader: skip the trailing incomplete frame of a live file."""
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise StructureError(
            "trajectory is missing or empty: {}".format(path))
    text = path.read_text()
    lines = text.splitlines()
    initial_frames, partial = frame_prefix(lines)
    if not initial_frames:
        raise StructureError(
            "{}: no complete frame found (first partial frame: {})".format(
                path, partial))
    frames = list(initial_frames)
    dropped = None
    while True:
        end = frames[-1][1]
        tmp = tempfile.NamedTemporaryFile(
            "w", suffix=".xyz", prefix="cs_coupling_", delete=False)
        try:
            tmp.write("\n".join(lines[:end]) + "\n")
            tmp.close()
            try:
                records, dedup = ub.load_trajectory(tmp.name, box)
            except ub.StrictError:
                if dropped is not None or len(frames) <= 1:
                    raise
                # The last complete-looking frame failed strict parsing: the
                # file was probably still being written when it was read.  Drop
                # exactly one frame, never mask deeper structural errors.
                dropped = {
                    "start_line": frames[-1][0] + 1,
                    "reason": "last frame failed strict parsing "
                              "(treated as still being written)",
                }
                frames = frames[:-1]
                continue
        finally:
            try:
                tmp.close()
            except Exception:
                pass
            try:
                Path(tmp.name).unlink()
            except OSError:
                pass
        break
    info = {
        "path": str(path),
        "sha256_of_read_text": __import__("hashlib").sha256(
            text.encode("utf-8")).hexdigest(),
        "raw_lines": len(lines),
        "complete_frames": len(records),
        "frame_candidates": len(initial_frames),
        "in_progress": partial is not None or dropped is not None,
        "skipped_partial_frame": partial,
        "dropped_strict_parse_frame": dropped,
        "note": ("trailing incomplete frames are skipped so a running CP2K job "
                 "can be inspected without waiting for the file to close; the "
                 "file itself is never modified"),
    }
    return records, dedup, info


def discover_from_dir(root):
    root = Path(root)
    if not root.is_dir():
        raise StructureError("--dir is not a directory: {}".format(root))
    trajectory = None
    for name in TRAJECTORY_CANDIDATES:
        if (root / name).is_file():
            trajectory = root / name
            break
    if trajectory is None:
        candidates = sorted(root.glob("*-pos-1.xyz"))
        if candidates:
            trajectory = candidates[0]
    if trajectory is None:
        raise StructureError(
            "no *-pos-1.xyz trajectory found in {}".format(root))
    energy = None
    for name in ENERGY_CANDIDATES:
        if (root / name).is_file():
            energy = root / name
            break
    if energy is None:
        candidates = sorted(root.glob("*.ener"))
        if candidates:
            energy = candidates[0]
    return trajectory, energy


def resolve_manifest_and_bias(root, explicit_manifest, explicit_biased):
    """Return (manifest dict or None, manifest path or None, biased, source)."""
    manifest_path = None
    if explicit_manifest:
        manifest_path = Path(explicit_manifest)
        if not manifest_path.is_file():
            raise StructureError("--manifest not found: {}".format(
                manifest_path))
    else:
        for name in MANIFEST_CANDIDATES:
            candidate = Path(root) / name
            if candidate.is_file():
                manifest_path = candidate
                break
    manifest = None
    biased = None
    source = None
    if manifest_path is not None:
        manifest = read_json(manifest_path, "input manifest")
        protocol = manifest.get("protocol")
        if protocol in PROTOCOL_BIASED:
            biased = True
            source = "manifest protocol {!r} ({})".format(
                protocol, manifest_path.name)
        else:
            checks = manifest.get("checks")
            output = manifest.get("output") if isinstance(
                manifest.get("output"), dict) else {}
            output_name = Path(str(output.get("path", ""))).name
            if (isinstance(checks, dict)
                    and checks.get("project") == "iface_unbiased") \
                    or output_name == "iface_unbiased.inp":
                biased = False
                source = "manifest identifies the unbiased input ({})".format(
                    manifest_path.name)
    if explicit_biased in ("true", "false"):
        biased = explicit_biased == "true"
        source = "CLI --biased {}".format(explicit_biased)
    elif explicit_biased not in ("auto",):
        raise StructureError(
            "--biased must be auto/true/false, got {!r}".format(
                explicit_biased))
    if biased is None:
        raise StructureError(
            "cannot determine the biased/unbiased provenance: pass --manifest "
            "or --biased true|false")
    return manifest, manifest_path, biased, source


def resolve_atom_indices(manifest, root):
    """Resolve 1-based C/S/O390/H393/water indices from atom_map when available."""
    atom_map = None
    source = None
    if isinstance(manifest, dict):
        seed_block = manifest.get("seed")
        if isinstance(seed_block, dict) and isinstance(
                seed_block.get("cluster_atom_map"), list):
            atom_map = seed_block["cluster_atom_map"]
            source = "input manifest seed.cluster_atom_map"
    if atom_map is None:
        for name in SEED_MANIFEST_CANDIDATES:
            candidate = Path(root) / name
            if candidate.is_file():
                seed_manifest = read_json(candidate, "seed manifest")
                cluster = seed_manifest.get("cluster")
                if isinstance(cluster, dict) and isinstance(
                        cluster.get("atom_map"), list):
                    atom_map = cluster["atom_map"]
                    source = "seed manifest cluster.atom_map ({})".format(name)
                    break
    indices = {
        "C": ub.INDEX_C + 1,
        "S": ub.INDEX_S + 1,
        "O_hydroxyl": ub.INDEX_O_HYDROXYL + 1,
        "H_hydroxyl": ub.INDEX_H_HYDROXYL + 1,
        "O_carbonyl": common.IDX_O_CARBONYL + 1,
        "water_O": None,
        "water_H1": None,
        "water_H2": None,
    }
    if atom_map is not None:
        by_role = {}
        for entry in atom_map:
            role = entry.get("role")
            source_1 = entry.get("source_index_1based")
            if role in ROLE_KEYS.values() and isinstance(source_1, int):
                by_role[role] = source_1
        for key, role in ROLE_KEYS.items():
            if role in by_role:
                indices[key] = by_role[role]
        source = (source or "atom map") + " (roles resolved)"
    return indices, source


def series(values):
    return ub.summarize([value for value in values if value is not None])


def longest_run(flags):
    best = 0
    current = 0
    best_end = None
    for index, flag in enumerate(flags):
        if flag:
            current += 1
            if current > best:
                best = current
                best_end = index
        else:
            current = 0
    return best, best_end


def first_true(records, flags):
    for record, flag in zip(records, flags):
        if flag:
            return {"step": record.step, "time_fs": record.time}
    return None


def trailing_runs(records, flags):
    """Trailing contiguous True segment -> (start index or None, count)."""
    index = len(records) - 1
    count = 0
    while index >= 0 and flags[index]:
        count += 1
        index -= 1
    start = None if count == 0 else index + 1
    return start, count


def analyze(trajectory_path, energy_path, box, manifest, manifest_path, biased,
            biased_source, indices, indices_source, relative_time=False):
    records, dedup_info, read_info = load_trajectory_in_progress(
        trajectory_path, box)
    if len(records[0].atoms) != ub.EXPECTED_ATOMS:
        raise StructureError(
            "{}: first frame has {} atoms, expected {}".format(
                trajectory_path, len(records[0].atoms), ub.EXPECTED_ATOMS))
    energy_records, energy_info, energy_error = ub.load_energy_optional(
        None if energy_path is None else Path(energy_path))
    if energy_records is None:
        temperature_lookup = lambda step, time_value: None  # noqa: E731
    else:
        temperature_lookup = ub.build_temperature_lookup(energy_records)
    records = [record._replace(
        temp=temperature_lookup(record.step, record.time))
        for record in records]
    time_origin = {"relative_time": False,
                   "note": "Absolute CP2K step/time are used."}
    if relative_time:
        records, time_origin = ub.normalize_time_origin(records)
        if energy_records is not None:
            energy_records = [record._replace(
                step=record.step - time_origin["step_offset"],
                time=record.time - time_origin["time_offset_fs"])
                for record in energy_records]

    sc_values = [record.sc for record in records]
    oh_values = [record.oh for record in records]
    cs_flags = [value < ub.S_C_FORMED_MAX_A for value in sc_values]
    oh_flags = [value > ub.O_H_BROKEN_MIN_A for value in oh_values]
    paper_flags = [cs and oh for cs, oh in zip(cs_flags, oh_flags)]

    first = records[0]
    last = records[-1]
    initial_unreacted = first.sc > ub.S_C_FORMED_MAX_A

    cs_trailing_index, cs_persistence = trailing_runs(records, cs_flags)
    cs_sustained = cs_persistence >= ub.MIN_PERSISTENCE_FRAMES
    cs_sustained_fs = (None if cs_trailing_index is None
                       else records[cs_trailing_index].time)
    cs_sustained_step = (None if cs_trailing_index is None
                         else records[cs_trailing_index].step)
    cs_within_3ps = (cs_sustained and cs_sustained_fs is not None
                     and cs_sustained_fs <= ub.MIN_COVERAGE_FS)
    cs_longest, cs_longest_end = longest_run(cs_flags)
    cs_first = first_true(records, cs_flags)

    paper_trailing_index, paper_persistence = trailing_runs(records,
                                                            paper_flags)
    paper_sustained = paper_persistence >= ub.MIN_PERSISTENCE_FRAMES
    paper_sustained_fs = (None if paper_trailing_index is None
                          else records[paper_trailing_index].time)
    paper_within_3ps = (paper_sustained and paper_sustained_fs is not None
                        and paper_sustained_fs <= ub.MIN_COVERAGE_FS)

    oh_longest, _oh_end = longest_run(oh_flags)
    oh_trailing_index, oh_persistence = trailing_runs(records, oh_flags)
    oh_sustained = oh_persistence >= ub.MIN_PERSISTENCE_FRAMES
    proton_transfer = {
        "occurred": any(oh_flags),
        "first_time_fs": (None if not any(oh_flags) else
                          first_true(records, oh_flags)["time_fs"]),
        "longest_run_frames": oh_longest,
        "sustained": oh_sustained,
        "sustained_fs": (None if oh_trailing_index is None
                         else records[oh_trailing_index].time),
        "persistence_frames": oh_persistence,
        "threshold_A": ub.O_H_BROKEN_MIN_A,
        "note": ("proton transfer is a SECONDARY observable for D18; it is "
                 "not part of the primary C-S criterion"),
    }

    water_distances = None
    if indices.get("water_O") is not None:
        d2 = []
        d3 = []
        d5 = []
        water_o = indices["water_O"] - 1
        water_h1 = indices["water_H1"] - 1
        o_carbonyl = indices["O_carbonyl"] - 1
        h_hydroxyl = indices["H_hydroxyl"] - 1
        for record in records:
            atoms = record.atoms
            d2.append(ub.periodic_distance(atoms[water_o][1:4],
                                           atoms[water_h1][1:4], box))
            d3.append(ub.periodic_distance(atoms[water_h1][1:4],
                                           atoms[o_carbonyl][1:4], box))
            d5.append(ub.periodic_distance(atoms[water_o][1:4],
                                           atoms[h_hydroxyl][1:4], box))
        water_distances = {
            "d2_water_OH_A": series(d2),
            "d3_waterH_O386_A": series(d3),
            "d5_waterO_H393_A": series(d5),
            "final": {"d2_A": d2[-1], "d3_A": d3[-1], "d5_A": d5[-1]},
            "min": {"d2_A": min(d2), "d3_A": min(d3), "d5_A": min(d5)},
            "note": ("d2/d3/d5 are the HMS water-relay diagnostics; they are "
                     "reported for interpretation only and are never used as "
                     "the primary D18 gate"),
        }

    coverage = ub.coverage_block(records)
    temperature = ub.temperature_block(
        records, energy_records, energy_error,
        cs_sustained_fs if cs_sustained else None)

    checks = {
        "initial_state": {
            "passed": initial_unreacted,
            "reason": ("first frame is unreacted (S-C {:.3f} A > {:.1f} A)"
                       .format(first.sc, ub.S_C_FORMED_MAX_A)
                       if initial_unreacted else
                       "first frame already has S-C {:.3f} A < {:.1f} A"
                       .format(first.sc, ub.S_C_FORMED_MAX_A)),
        },
        "trajectory_coverage": coverage,
        "cs_sustained": {
            "passed": cs_sustained,
            "reason": (
                "trailing C-S bonded segment starts at {:.3f} fs and covers "
                "{} frames".format(cs_sustained_fs, cs_persistence)
                if cs_sustained else
                "no trailing C-S bonded segment of at least {} frames "
                "(longest run {} frames)".format(
                    ub.MIN_PERSISTENCE_FRAMES, cs_longest)),
        },
        "cs_within_3ps": {
            "passed": cs_within_3ps,
            "reason": (
                "sustained C-S segment starts at {:.3f} fs, within {:.0f} fs"
                .format(cs_sustained_fs, ub.MIN_COVERAGE_FS)
                if cs_within_3ps else
                ("no sustained C-S segment" if cs_sustained_fs is None else
                 "sustained C-S segment starts at {:.3f} fs, after {:.0f} fs"
                 .format(cs_sustained_fs, ub.MIN_COVERAGE_FS))),
        },
        "temperature": temperature,
    }

    if not coverage["passed"]:
        verdict, exit_code = "CS_COUPLING_COVERAGE_SHORT", 5
    elif not initial_unreacted:
        verdict, exit_code = "CS_COUPLING_INITIAL_STATE_NOT_UNREACTED", 4
    elif not cs_sustained:
        verdict, exit_code = VERDICTS[4], 4
    elif not temperature["passed"]:
        verdict, exit_code = "CS_COUPLING_SUSTAINED_TEMPERATURE_FAIL", 4
    else:
        verdict, exit_code = VERDICTS[0], 0
    failed_checks = [name for name in checks if not checks[name]["passed"]]
    summary = {
        "status": "ok",
        "protocol": PROTOCOL,
        "criterion": CRITERION,
        "criterion_is_paper_comparable": False,
        "not_paper_criterion": True,
        "criterion_note": (
            "D18 primary criterion = C-S bond only (S-C < {:.1f} A sustained "
            ">= {} frames).  This is the user-specified primary event on "
            "2026-10-28; it is NOT the paper-comparable criterion.  The "
            "original S-C AND proton-transfer criterion is reported separately "
            "under paper_criterion and must be quoted alongside.".format(
                ub.S_C_FORMED_MAX_A, ub.MIN_PERSISTENCE_FRAMES)),
        "biased": bool(biased),
        "biased_source": biased_source,
        "input_manifest": (None if manifest_path is None else {
            "path": str(manifest_path),
            "sha256": sha256(manifest_path),
            "protocol": (manifest or {}).get("protocol"),
        }),
        "source_trajectory": str(trajectory_path),
        "source_energy": (None if energy_path is None else str(energy_path)),
        "box_angstrom": list(box),
        "atom_count": ub.EXPECTED_ATOMS,
        "atom_indices": dict(indices, source=indices_source),
        "trajectory_read": read_info,
        "deduplication": dedup_info,
        "energy": energy_info,
        "time_origin": time_origin,
        "criteria": {
            "s_c_formed_max_A": ub.S_C_FORMED_MAX_A,
            "o_h_broken_min_A": ub.O_H_BROKEN_MIN_A,
            "min_persistence_frames": ub.MIN_PERSISTENCE_FRAMES,
            "min_coverage_fs": ub.MIN_COVERAGE_FS,
            "max_adjacent_gap_fs": ub.MAX_ADJACENT_GAP_FS,
            "min_unique_frames": ub.MIN_UNIQUE_FRAMES,
            "source": ("thresholds and persistence rule imported from "
                       "analyze_unbiased (not re-typed)"),
        },
        "primary": {
            "criterion": CRITERION,
            "s_c_min_A": min(sc_values),
            "s_c_final_A": last.sc,
            "s_c_first_below_2A": cs_first,
            "frames_below_2A": sum(1 for flag in cs_flags if flag),
            "longest_run_frames": cs_longest,
            "longest_run_end_step": (None if cs_longest_end is None
                                     else records[cs_longest_end].step),
            "sustained": cs_sustained,
            "sustained_fs": cs_sustained_fs,
            "sustained_step": cs_sustained_step,
            "sustained_duration_fs": (None if cs_trailing_index is None
                                      else last.time - cs_sustained_fs),
            "persistence_frames": cs_persistence,
            "minimum_persistence_frames": ub.MIN_PERSISTENCE_FRAMES,
            "sustained_within_3ps": cs_within_3ps,
            "final_frame": {"step": last.step, "time_fs": last.time,
                            "s_c_A": last.sc, "o_h_A": last.oh},
        },
        "paper_criterion": {
            "criterion": "paper_comparable_s_c_and_proton_transfer",
            "definition": ("S-C < {:.1f} A AND O-H > {:.1f} A, same "
                           "trailing-{}-frame persistence rule".format(
                               ub.S_C_FORMED_MAX_A, ub.O_H_BROKEN_MIN_A,
                               ub.MIN_PERSISTENCE_FRAMES)),
            "sustained": paper_sustained,
            "sustained_fs": paper_sustained_fs,
            "persistence_frames": paper_persistence,
            "sustained_within_3ps": paper_within_3ps,
            "note": ("reported for comparison with the paper; it is NOT the "
                     "primary D18 criterion and did not gate this verdict"),
        },
        "secondary": {
            "o_h": series(oh_values),
            "o_h_final_A": last.oh,
            "proton_transfer": proton_transfer,
            "water_distances": water_distances,
        },
        "checks": checks,
        "failed_checks": failed_checks,
        "failure_reasons": [checks[name]["reason"] for name in failed_checks],
        "verdict": verdict,
        "exit_code": exit_code,
        "semantics": {
            "question_answered": "did a C-S bond form and persist?",
            "question_not_answered": [
                "did proton transfer happen (reported as secondary)",
                "is this the paper-comparable product definition",
                "was the reaction spontaneous (biased runs cannot answer this)",
            ],
            "reporting_rule": (
                "quote the D18 primary C-S result together with the "
                "paper_criterion result and the biased flag; guided runs are "
                "mechanism demonstrations only"),
        },
    }
    return summary, records, exit_code


def write_outputs(summary, records, json_path, tsv_path, d2, d3, d5):
    json_path = Path(json_path)
    tsv_path = Path(tsv_path)
    if json_path.parent != Path(""):
        json_path.parent.mkdir(parents=True, exist_ok=True)
    if tsv_path.parent != Path(""):
        tsv_path.parent.mkdir(parents=True, exist_ok=True)
    json_text = json.dumps(summary, indent=2) + "\n"
    std_line = (
        "CS_COUPLING_VERDICT verdict={} criterion={} biased={} "
        "s_c_min_A={:.4f} sustained_frames={} paper_criterion_sustained={} "
        "exit_code={}".format(
            summary["verdict"], summary["criterion"], summary["biased"],
            summary["primary"]["s_c_min_A"],
            summary["primary"]["persistence_frames"],
            summary["paper_criterion"]["sustained"],
            summary["exit_code"]))
    require_cs_semantics(json_text + std_line, summary)
    json_path.write_text(json_text)
    with tsv_path.open("w") as handle:
        handle.write("step\ttime_fs\ts_c_A\to_h_A\td2_A\td3_A\td5_A\ttemp_K\t"
                     "cs_bond_frame\tpaper_product_frame\n")
        for index, record in enumerate(records):
            temperature = ("" if record.temp is None
                           else "{:.3f}".format(record.temp))
            handle.write("{}\t{:.3f}\t{:.5f}\t{:.5f}\t{}\t{}\t{}\t{}\t{}\t{}\n"
                         .format(record.step, record.time, record.sc,
                                 record.oh,
                                 "" if d2 is None else "{:.5f}".format(d2[index]),
                                 "" if d3 is None else "{:.5f}".format(d3[index]),
                                 "" if d5 is None else "{:.5f}".format(d5[index]),
                                 temperature,
                                 1 if record.sc < ub.S_C_FORMED_MAX_A else 0,
                                 1 if record.product else 0))
    return std_line


def require_cs_semantics(text, summary):
    for token in FORBIDDEN_OUTPUT_TOKENS:
        if token in text:
            raise StructureError(
                "D18 output contains the paper-criterion verdict token {!r}; "
                "refusing to mix the two semantics".format(token))
    if summary.get("criterion") != CRITERION:
        raise StructureError(
            "D18 output must declare criterion={!r}".format(CRITERION))
    if summary.get("criterion_is_paper_comparable") is not False:
        raise StructureError(
            "D18 output must record criterion_is_paper_comparable=false")
    if not isinstance(summary.get("paper_criterion"), dict):
        raise StructureError(
            "D18 output must also report the paper-comparable criterion")
    if not isinstance(summary.get("biased"), bool):
        raise StructureError("D18 output must carry a boolean biased field")
    if not str(summary.get("verdict", "")).startswith("CS_COUPLING_"):
        raise StructureError(
            "D18 verdict must be a CS_COUPLING_* string, got {!r}".format(
                summary.get("verdict")))


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=("D18 C-S coupling analysis: primary criterion = C-S bond "
                     "only (S-C < 2.0 A sustained); the paper-comparable "
                     "S-C+proton-transfer criterion is reported separately."),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("trajectory", nargs="?",
                        help="XYZ trajectory (i=<step> time=<value> comments); "
                             "omit when using --dir")
    parser.add_argument("--dir", default=None,
                        help="BOMD directory; auto-discovers trajectory/.ener/"
                             "input manifest")
    parser.add_argument("--energy", default=None,
                        help="CP2K .ener file for the temperature evidence")
    parser.add_argument("--manifest", default=None,
                        help="input manifest (IFACE_STRONG_CS_INPUT.json / "
                             "IFACE_STEERED_INPUT.json / IFACE_UNBIASED_INPUT."
                             "json); auto-discovers when omitted")
    parser.add_argument("--biased", default="auto",
                        choices=("auto", "true", "false"),
                        help="biased provenance override (default: read from "
                             "the manifest)")
    parser.add_argument("--box", type=ub.box_argument, default=ub.DEFAULT_BOX,
                        help="orthorhombic cell lengths in A, default 17,17,47")
    parser.add_argument("--json", default=None,
                        help="output JSON (default: CS_COUPLING.json next to "
                             "the trajectory / in --dir)")
    parser.add_argument("--tsv", default=None,
                        help="output TSV (default: CS_COUPLING.tsv next to "
                             "the trajectory / in --dir)")
    parser.add_argument("--relative-time", action="store_true",
                        help="shift step/time to the first analyzed frame for "
                             "a segmented run; raw files and thresholds are "
                             "unchanged")
    args = parser.parse_args(argv)

    try:
        if args.dir and args.trajectory:
            raise StructureError("pass either a trajectory or --dir, not both")
        if args.dir:
            trajectory_path, discovered_energy = discover_from_dir(args.dir)
            root = Path(args.dir)
            energy_path = (Path(args.energy) if args.energy
                           else discovered_energy)
        else:
            if not args.trajectory:
                raise StructureError("pass a trajectory path or --dir")
            trajectory_path = Path(args.trajectory)
            root = trajectory_path.parent
            energy_path = Path(args.energy) if args.energy else None
        json_path = (Path(args.json) if args.json is not None
                     else root / DEFAULT_JSON)
        tsv_path = (Path(args.tsv) if args.tsv is not None
                    else root / DEFAULT_TSV)
        manifest, manifest_path, biased, biased_source = \
            resolve_manifest_and_bias(root, args.manifest, args.biased)
        indices, indices_source = resolve_atom_indices(manifest, root)
        summary, records, exit_code = analyze(
            trajectory_path, energy_path, args.box, manifest, manifest_path,
            biased, biased_source, indices, indices_source,
            relative_time=args.relative_time)
        water = summary["secondary"]["water_distances"]
        if water is None:
            d2 = d3 = d5 = None
        else:
            d2, d3, d5 = [], [], []
            water_o = indices["water_O"] - 1
            water_h1 = indices["water_H1"] - 1
            o_carbonyl = indices["O_carbonyl"] - 1
            h_hydroxyl = indices["H_hydroxyl"] - 1
            for record in records:
                atoms = record.atoms
                d2.append(ub.periodic_distance(atoms[water_o][1:4],
                                               atoms[water_h1][1:4], args.box))
                d3.append(ub.periodic_distance(atoms[water_h1][1:4],
                                               atoms[o_carbonyl][1:4],
                                               args.box))
                d5.append(ub.periodic_distance(atoms[water_o][1:4],
                                               atoms[h_hydroxyl][1:4],
                                               args.box))
        std_line = write_outputs(summary, records, json_path, tsv_path,
                                 d2, d3, d5)
    except (StructureError, RuntimeError, OSError, ValueError) as exc:
        sys.stderr.write("analyze_cs_coupling: error: {}\n".format(exc))
        return 2
    sys.stdout.write(std_line + "\n")
    sys.stdout.write("CS_COUPLING_JSON " + str(json_path) + "\n")
    return exit_code


if __name__ == "__main__":
    sys.exit(main())

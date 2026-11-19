#!/usr/bin/env python3
"""Factual per-frame analysis of an unbiased Li 2025 event replica.

Reads a CP2K XYZ trajectory (296 atoms, the fixed Li 2025 order) and reports,
for every frame: d(S-C), d(C=O), both d(C-H), the three d(S-O), and the nearest
water hydrogen to the carbonyl oxygen together with that water's covalent O-H.

Summary facts (no verdict):
  * first/last/min/max of each distance,
  * frames with d(S-C) < 2.0 A and the longest contiguous such run,
  * first sustained dissociation: d(S-C) > 2.5 A for 10 consecutive frames
    (a project recording threshold, not a literature criterion),
  * temperature coverage from the CP2K .ener file, when supplied.

This tool deliberately emits no PASS/FAIL verdict: the output says what the
distances did, and states explicitly that this is unbiased-replica stability
evidence, not a free energy and not proof of spontaneous formation.

Usage:
  python3 analyze_event_unbiased.py --trajectory run-pos-1.xyz \
      --energy run-1.ener --label event_unbiased_s1 --box 14.2296,14.2296,14.2296 \
      --out EVENT_UNBIASED_s1.json
"""
import argparse
import hashlib
import json
import math
import re
from pathlib import Path

ATOM_COUNT = 296
EXPECTED_ORDER = ["H", "H", "O"] * 96 + ["C", "O", "H", "H", "S", "O", "O", "O"]
C_INDEX = 288
O_CARBONYL = 289
H_FORMYL = (290, 291)
S_INDEX = 292
O_SULFONATE = (293, 294, 295)
WATER_H = list(range(0, 288, 3)) + list(range(1, 288, 3))


def sha256(path):
    digest = hashlib.sha256()
    with open(str(path), "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def distance(a, b, box):
    delta = [a[i] - b[i] for i in range(3)]
    delta = [delta[i] - box[i] * round(delta[i] / box[i]) for i in range(3)]
    return math.sqrt(sum(value * value for value in delta))


def iter_frames(path):
    comment_pattern = re.compile(r"i\s*=\s*(\d+)[,\s]+time\s*=\s*([-0-9.]+)")
    with open(str(path), "r") as handle:
        while True:
            count_line = handle.readline()
            if not count_line:
                return
            if not count_line.strip():
                continue
            try:
                atom_count = int(count_line.strip())
            except ValueError:
                raise SystemExit("bad atom-count line in %s: %r" % (path, count_line))
            if atom_count != ATOM_COUNT:
                raise SystemExit("%s has %d atoms, expected %d" % (path, atom_count, ATOM_COUNT))
            comment = handle.readline().rstrip("\n")
            elements = []
            coords = []
            for _ in range(atom_count):
                parts = handle.readline().split()
                if len(parts) != 4:
                    raise SystemExit("bad atom line in %s: %r" % (path, parts))
                elements.append(parts[0])
                coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
            if elements != EXPECTED_ORDER:
                raise SystemExit("atom order in %s is not the Li 2025 order" % path)
            match = comment_pattern.search(comment)
            if not match:
                raise SystemExit("frame comment lacks i/time in %s: %r" % (path, comment))
            yield int(match.group(1)), float(match.group(2)), coords


def frame_metrics(coords, box):
    carbon = coords[C_INDEX]
    oxygen = coords[O_CARBONYL]
    sulfur = coords[S_INDEX]
    so = [distance(sulfur, coords[index], box) for index in O_SULFONATE]

    best = (1.0e9, None, None)
    for hydrogen in WATER_H:
        value = distance(oxygen, coords[hydrogen], box)
        if value < best[0]:
            best = (value, hydrogen, None)
    # Owning water oxygen: the O of the same H,H,O group.
    hydrogen = best[1]
    water_oxygen = hydrogen - (hydrogen % 3) + 2
    covalent_oh = distance(coords[water_oxygen], coords[hydrogen], box)

    return {
        "sc": distance(sulfur, carbon, box),
        "co": distance(carbon, oxygen, box),
        "ch1": distance(carbon, coords[H_FORMYL[0]], box),
        "ch2": distance(carbon, coords[H_FORMYL[1]], box),
        "so1": so[0],
        "so2": so[1],
        "so3": so[2],
        "so_mean": sum(so) / 3.0,
        "nearest_water_h_1based": hydrogen + 1,
        "nearest_water_h_distance": best[0],
        "water_o_1based": water_oxygen + 1,
        "water_oh_distance": covalent_oh,
    }


def scalar_summary(values):
    return {
        "first": values[0],
        "last": values[-1],
        "min": min(values),
        "max": max(values),
        "mean": sum(values) / float(len(values)),
    }


def longest_run(flags):
    best_length = 0
    best_start = None
    best_end = None
    start = None
    for index, flag in enumerate(flags):
        if flag and start is None:
            start = index
        if not flag and start is not None:
            length = index - start
            if length > best_length:
                best_length, best_start, best_end = length, start, index - 1
            start = None
    if start is not None:
        length = len(flags) - start
        if length > best_length:
            best_length, best_start, best_end = length, start, len(flags) - 1
    return best_length, best_start, best_end


def first_sustained(values, threshold, run_length):
    for index in range(0, len(values) - run_length + 1):
        if all(value > threshold for value in values[index:index + run_length]):
            return index
    return None


def temperature_summary(path, time_range=None):
    """CP2K .ener columns are Step, Time, Kin, Temp, Pot, Cons, Used (7 fields)."""
    if path is None:
        return None
    path = Path(path)
    if not path.is_file():
        return None
    header_seen = False
    temps = []
    with open(str(path), "r") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                if "Temp[K]" in stripped:
                    header_seen = True
                continue
            parts = stripped.split()
            if not header_seen or len(parts) != 7:
                continue
            try:
                step = int(parts[0])
                time_fs = float(parts[1])
                temperature = float(parts[3])
            except ValueError:
                continue
            if time_range is not None and not (time_range[0] <= time_fs <= time_range[1]):
                continue
            temps.append((step, time_fs, temperature))
    if not temps:
        return {"source": str(path), "rows": 0}
    target = 298.15
    values = [item[2] for item in temps]
    within = sum(1 for value in values if abs(value - target) <= 25.0)
    return {
        "source": str(path),
        "rows": len(values),
        "first_step": temps[0][0],
        "last_step": temps[-1][0],
        "first_k": values[0],
        "last_k": values[-1],
        "min_k": min(values),
        "max_k": max(values),
        "mean_k": sum(values) / float(len(values)),
        "within_25k_of_298.15_fraction": within / float(len(values)),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trajectory", required=True)
    parser.add_argument("--energy", default=None)
    parser.add_argument("--out", required=True)
    parser.add_argument("--label", default=None)
    parser.add_argument("--box", required=True,
                        help="cell lengths in A, comma separated: 14.2296,14.2296,14.2296 "
                             "for the bulk cube or 15,15,30 for the interface cell. Required "
                             "on purpose: a default would silently apply the bulk cell to an "
                             "interface trajectory and corrupt every minimum-image distance.")
    parser.add_argument("--bond-threshold", type=float, default=2.0)
    parser.add_argument("--dissociation-threshold", type=float, default=2.5)
    parser.add_argument("--dissociation-frames", type=int, default=10)
    args = parser.parse_args()

    box = [float(value) for value in args.box.split(",")]
    if len(box) != 3:
        raise SystemExit("--box needs three comma-separated lengths")

    trajectory = Path(args.trajectory)
    rows = []
    for step, time_fs, coords in iter_frames(trajectory):
        metrics = frame_metrics(coords, box)
        metrics["step"] = step
        metrics["time_fs"] = time_fs
        rows.append(metrics)
    if not rows:
        raise SystemExit("no frames read from %s" % trajectory)

    sc = [row["sc"] for row in rows]
    flags = [value < args.bond_threshold for value in sc]
    run_length, run_start, run_end = longest_run(flags)
    dissociation_index = first_sustained(sc, args.dissociation_threshold,
                                         args.dissociation_frames)

    summary = {
        "sc": scalar_summary(sc),
        "co": scalar_summary([row["co"] for row in rows]),
        "ch1": scalar_summary([row["ch1"] for row in rows]),
        "ch2": scalar_summary([row["ch2"] for row in rows]),
        "so_mean": scalar_summary([row["so_mean"] for row in rows]),
        "nearest_water_h_distance": scalar_summary([row["nearest_water_h_distance"]
                                                    for row in rows]),
        "water_oh_distance": scalar_summary([row["water_oh_distance"] for row in rows]),
    }
    summary["sc"]["frames_below_threshold"] = sum(1 for flag in flags if flag)
    summary["sc"]["bond_threshold_angstrom"] = args.bond_threshold
    summary["sc"]["longest_below_threshold_frames"] = run_length
    summary["sc"]["longest_below_threshold_window"] = (
        None if run_start is None else
        {"start_index": run_start, "end_index": run_end,
         "start_step": rows[run_start]["step"], "end_step": rows[run_end]["step"],
         "start_time_fs": rows[run_start]["time_fs"],
         "end_time_fs": rows[run_end]["time_fs"]})
    if dissociation_index is None:
        summary["sc"]["first_sustained_dissociation"] = None
    else:
        end = dissociation_index + args.dissociation_frames - 1
        summary["sc"]["first_sustained_dissociation"] = {
            "criterion": "d(S-C) > %.2f A for %d consecutive frames"
                         % (args.dissociation_threshold, args.dissociation_frames),
            "index": dissociation_index,
            "step": rows[dissociation_index]["step"],
            "time_fs": rows[dissociation_index]["time_fs"],
            "window_end_step": rows[end]["step"],
            "window_end_time_fs": rows[end]["time_fs"],
        }

    result = {
        "label": args.label or trajectory.stem,
        "trajectory": str(trajectory),
        "trajectory_sha256": sha256(trajectory),
        "atom_count": ATOM_COUNT,
        "box_angstrom": box,
        "frames": {
            "count": len(rows),
            "first_step": rows[0]["step"],
            "last_step": rows[-1]["step"],
            "first_time_fs": rows[0]["time_fs"],
            "last_time_fs": rows[-1]["time_fs"],
        },
        "summary": summary,
        "temperature": temperature_summary(args.energy,
                                           (rows[0]["time_fs"], rows[-1]["time_fs"])),
        "per_frame": rows,
        "interpretation": {
            "what_this_is": ("unbiased BOMD replica stability evidence: distance facts "
                             "after removing metadynamics"),
            "what_this_is_not": [
                "not a free energy, PMF, barrier or rate",
                "not proof of spontaneous formation (the start frame came from a metaD run)",
                "not a Li three-replica statistic",
                "not an electronic bonding analysis",
            ],
            "verdict": None,
            "note": ("no pass/fail verdict is asserted; the project's C-S bonding gate "
                     "requires the three-layer review in AGENTS.md section 7"),
        },
    }

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

    tsv_path = out_path.with_name(out_path.stem + "_FRAMES.tsv")
    keys = ["step", "time_fs", "sc", "co", "ch1", "ch2", "so1", "so2", "so3",
            "so_mean", "nearest_water_h_1based", "nearest_water_h_distance",
            "water_o_1based", "water_oh_distance"]
    with open(str(tsv_path), "w") as handle:
        handle.write("\t".join(keys) + "\n")
        for row in rows:
            handle.write("\t".join("%.6f" % row[key] if isinstance(row[key], float)
                                   else str(row[key]) for key in keys) + "\n")

    print(json.dumps({
        "label": result["label"],
        "frames": result["frames"],
        "sc": summary["sc"],
        "temperature": result["temperature"],
        "out": str(out_path),
        "tsv": str(tsv_path),
    }, indent=2))


if __name__ == "__main__":
    main()

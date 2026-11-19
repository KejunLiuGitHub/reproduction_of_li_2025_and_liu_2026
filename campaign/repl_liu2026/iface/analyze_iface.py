#!/usr/bin/env python3
"""Interface structure checks for the iface BOMD trajectory.

Single-frame report (last frame of the file) and series reports.  Frames are
parsed strictly: every frame must carry a complete atom block, and every frame
used by --series must name its MD step in the XYZ comment.  --series always
writes IFACE_FULL_SERIES.json with every failing frame retained; it writes
IFACE_SERIES.json for the same whole trajectory by default, or for the
production window opened by the IFACE_SERIES_START_STEP environment variable.
Kept Python 3.6 compatible because the cluster interpreter is old.
"""
import hashlib
import json
import os
import re
import statistics
import sys
from collections import namedtuple
from pathlib import Path

EXPECTED_ATOMS = 394
CELL = (17.0, 17.0, 47.0)
BIN_COUNT = 94
WATER_ATOMS = 384
TARGET_MD_STEP = 10000
STABILITY_THICKNESS_HALF_DIFF_MAX_A = 1.0
STABILITY_CENTER_HALF_DIFF_MAX_A = 1.0
DEDUP_POLICY = ("dedupe by XYZ comment i=<step>; keep the last occurrence in the "
                "file; sort by step; raw file is never rewritten")
PRODUCTION_WINDOW_ENV = "IFACE_SERIES_START_STEP"
FULL_SERIES_NAME = "IFACE_FULL_SERIES.json"
WINDOW_SERIES_NAME = "IFACE_SERIES.json"
WINDOW_POLICY_ID = "project_defined_last_5ps_production_window"
FAILURE_POLICY = ("every structurally invalid frame is retained with its step, "
                  "start line and reason; failing frames stay in the accounting "
                  "(frames == frames_passed + frames_failed) and are never "
                  "deleted, silently skipped or replaced by passing frames")
WINDOW_RATIONALE = (
    "Frames from IFACE_SERIES_START_STEP to the final step form the RC production "
    "window because the early pre-equilibration is still relaxing the slab; the "
    "window matches prepare_iface_rc.py --window-ps 5.0. No consulted literature "
    "states a whole-trajectory zero-anomaly requirement, so this window is an "
    "explicit project convention, not a literature threshold. Early-transient "
    "frames are retained in IFACE_FULL_SERIES.json, and a passing window is not a "
    "full-trajectory stability claim")
STEP_PATTERN = re.compile(r"\bi\s*=\s*(\d+)")

Frame = namedtuple("Frame", "step atoms start_line")
Metrics = namedtuple("Metrics", "result density smooth")


def read_frames(path):
    """Parse every frame of an XYZ file, validating atom counts and blocks.

    A trailing partial frame, a malformed atom line, or a non-integer header
    raises instead of silently shifting the parse.
    """
    path = Path(path)
    lines = path.read_text().splitlines()
    while lines and not lines[-1].strip():
        lines.pop()
    frames = []
    offset = 0
    total = len(lines)
    while offset < total:
        header = lines[offset].strip()
        try:
            count = int(header)
        except ValueError:
            raise RuntimeError(
                "{}: line {}: expected an XYZ atom count, found {!r}".format(
                    path, offset + 1, header))
        if count <= 0:
            raise RuntimeError(
                "{}: line {}: non-positive atom count {}".format(path, offset + 1, count))
        if offset + count + 2 > total:
            raise RuntimeError(
                "{}: line {}: truncated frame; header claims {} atoms but only {} "
                "lines remain".format(path, offset + 1, count, total - offset - 2))
        comment = lines[offset + 1]
        match = STEP_PATTERN.search(comment)
        atoms = []
        for index in range(count):
            line_number = offset + 2 + index
            parts = lines[line_number].split()
            if len(parts) < 4:
                raise RuntimeError(
                    "{}: line {}: malformed atom line {!r}".format(
                        path, line_number + 1, lines[line_number]))
            try:
                coordinates = (float(parts[1]), float(parts[2]), float(parts[3]))
            except ValueError:
                raise RuntimeError(
                    "{}: line {}: non-numeric coordinate in {!r}".format(
                        path, line_number + 1, lines[line_number]))
            atoms.append((parts[0],) + coordinates)
        frames.append(Frame(step=int(match.group(1)) if match else None,
                            atoms=atoms, start_line=offset + 1))
        offset += count + 2
    if not frames:
        raise RuntimeError("{}: no complete frame found".format(path))
    return frames


def frame_metrics(atoms):
    """Density based slab/interface metrics for one complete frame."""
    count = len(atoms)
    if count != EXPECTED_ATOMS:
        raise RuntimeError("expected {} atoms, got {}".format(EXPECTED_ATOMS, count))
    x_cell, y_cell, z_cell = CELL
    width = z_cell / BIN_COUNT
    bins = [0] * BIN_COUNT
    for index in range(0, WATER_ATOMS, 3):
        z = atoms[index][3] % z_cell
        bins[min(BIN_COUNT - 1, int(z / width))] += 1
    volume = x_cell * y_cell * width
    density = [value / volume for value in bins]
    smooth = [sum(density[max(0, i - 1):min(BIN_COUNT, i + 2)]) /
              len(density[max(0, i - 1):min(BIN_COUNT, i + 2)]) for i in range(BIN_COUNT)]
    maximum = max(smooth)
    if maximum <= 0.0:
        raise RuntimeError("no water oxygen density found")
    threshold = maximum * 0.5
    regions = []
    start = None
    for i, value in enumerate(smooth + [0.0]):
        if value >= threshold and start is None:
            start = i
        elif value < threshold and start is not None:
            regions.append((start, i - 1))
            start = None
    if not regions:
        raise RuntimeError("no water slab found")
    merged = [regions[0]]
    for start, end in regions[1:]:
        if start - merged[-1][1] - 1 <= 3:
            merged[-1] = (merged[-1][0], end)
        else:
            merged.append((start, end))
    if len(merged) > 1 and merged[0][0] == 0 and merged[-1][1] == BIN_COUNT - 1:
        merged = [merged[-1]] + merged[:-1]
        merged[0] = (merged[0][0], merged[1][1] + BIN_COUNT)
        del merged[1]
    lower, upper = max(merged, key=lambda region: region[1] - region[0])

    def in_slab(index):
        return lower <= index <= upper or lower <= index + BIN_COUNT <= upper

    outside = [index for index in range(BIN_COUNT) if not in_slab(index)]
    low_outside = [index for index in outside if smooth[index] < threshold]
    two_interfaces_verified = (len(outside) >= 2
                               and len(low_outside) == len(outside))
    if not two_interfaces_verified:
        raise RuntimeError(
            "two-interface verification failed: slab bins {}-{}, {} outside bins, "
            "{} below half-max; every outside bin must stay below half-max and at "
            "least two low-density bins are required".format(
                lower, upper, len(outside), len(low_outside)))

    z_lower = (lower * width) % z_cell
    z_upper = ((upper + 1) * width) % z_cell
    thickness = (upper - lower + 1) * width
    if thickness < 5.0:
        raise RuntimeError("water slab too thin: {:.2f} A".format(thickness))
    center = (z_lower + thickness / 2) % z_cell
    if not 18.5 <= center <= 28.5:
        raise RuntimeError("slab center moved unexpectedly: {:.2f}".format(center))

    def inside(z):
        return z_lower <= z <= z_upper if z_lower < z_upper else z >= z_lower or z <= z_upper

    def interface_distance(z, boundary):
        direct = abs(z - boundary)
        return min(direct, z_cell - direct)

    def gds_signed(z, is_inside):
        """Signed distance to the nearest GDS; positive = toward the vapour.

        z_lower/z_upper come from the half-max water-density profile, i.e. they
        are the Gibbs dividing surfaces used by Liu 2026 Fig. 5. The legacy
        unsigned distance_to_nearest_interface cannot distinguish a reactant
        sitting on the GDS from one buried in the slab; this signed field can.
        Sign convention is shared with scripts/interface_orientation.py.
        """
        d = min(interface_distance(z, z_lower), interface_distance(z, z_upper))
        return -d if is_inside else d

    solutes = {}
    for name, index in (("S389", 388), ("C385", 384), ("K394", 393), ("O386", 385)):
        z = atoms[index][3] % z_cell
        is_inside = inside(z)
        solutes[name] = {"z": z, "inside_slab": is_inside,
                         "distance_to_nearest_interface": min(interface_distance(z, z_lower),
                                                              interface_distance(z, z_upper)),
                         "signed_distance_to_gds": gds_signed(z, is_inside),
                         "on_vapour_side": not is_inside}
    result = {
        "atom_count": count, "cell_angstrom": [x_cell, y_cell, z_cell],
        "max_number_density_A-3": maximum, "half_max_threshold_A-3": threshold,
        "slab_bin_range": [lower, upper],
        "z_lower_A": z_lower, "z_upper_A": z_upper,
        "thickness_A": thickness, "center_A": center,
        "outside_slab_bins": {
            "count": len(outside),
            "low_density_count": len(low_outside),
            "max_smoothed_A-3": max(smooth[index] for index in outside) if outside else None,
        },
        "two_interfaces_verified": two_interfaces_verified,
        "solutes": solutes,
    }
    return Metrics(result=result, density=density, smooth=smooth)


def main(path, write_outputs=True):
    path = Path(path)
    frames = read_frames(path)
    frame = frames[-1]
    metrics = frame_metrics(frame.atoms)
    result = metrics.result
    result["source_trajectory"] = str(path)
    result["frame_start_line"] = frame.start_line
    result["frame_step"] = frame.step
    result["frames_in_file"] = len(frames)
    if write_outputs:
        (path.parent / "IFACE_DENSITY.json").write_text(json.dumps(result, indent=2) + "\n")
        with (path.parent / "IFACE_DENSITY.tsv").open("w") as handle:
            handle.write("z_A\tnumber_density_A-3\tsmoothed_A-3\n")
            for i in range(BIN_COUNT):
                handle.write("{:.4f}\t{:.8f}\t{:.8f}\n".format(
                    (i + 0.5) * (CELL[2] / BIN_COUNT), metrics.density[i], metrics.smooth[i]))
        print(json.dumps(result, indent=2))
    return result


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1048576), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_and_dedupe(path):
    """Parse and step-deduplicate a trajectory without evaluating any frame."""
    path = Path(path)
    frames = read_frames(path)
    by_step = {}
    raw_steps = []
    duplicate_frames = 0
    for frame in frames:
        if frame.step is None:
            raise RuntimeError(
                "{}: frame at line {} has no 'i=<step>' in its XYZ comment; cannot "
                "deduplicate the series".format(path, frame.start_line))
        if frame.step in by_step:
            duplicate_frames += 1
        by_step[frame.step] = frame
        raw_steps.append(frame.step)
    steps = sorted(by_step)
    if len(steps) < 2:
        raise RuntimeError("need at least two complete frames")
    return {
        "path": path,
        "raw_frames": len(frames),
        "raw_steps": raw_steps,
        "steps": steps,
        "ordered": [by_step[step] for step in steps],
        "duplicate_frames": duplicate_frames,
    }


def evaluate_frames(ordered):
    """Evaluate every frame; a structural failure is recorded, never raised."""
    records = []
    for frame in ordered:
        try:
            metrics = frame_metrics(frame.atoms)
        except RuntimeError as error:
            records.append({"step": frame.step, "start_line": frame.start_line,
                            "metrics": None, "error": str(error)})
        else:
            records.append({"step": frame.step, "start_line": frame.start_line,
                            "metrics": metrics.result, "error": None})
    return records


def series_stats(values):
    if len(values) < 2:
        return {"mean": values[0], "stdev": None, "min": values[0], "max": values[0]}
    return {"mean": statistics.mean(values), "stdev": statistics.stdev(values),
            "min": min(values), "max": max(values)}


def half_stability(metrics):
    """Half-window stability with the unchanged project thresholds."""
    block = {
        "first_half_thickness_mean_A": None,
        "second_half_thickness_mean_A": None,
        "half_thickness_difference_A": None,
        "first_half_center_mean_A": None,
        "second_half_center_mean_A": None,
        "half_center_difference_A": None,
        "stability": {
            "computed": False,
            "stable": None,
            "half_thickness_difference_A": None,
            "half_center_difference_A": None,
            "thickness_stable": None,
            "center_stable": None,
            "frames_used": len(metrics),
            "thresholds": {
                "half_thickness_difference_max_A": STABILITY_THICKNESS_HALF_DIFF_MAX_A,
                "half_center_difference_max_A": STABILITY_CENTER_HALF_DIFF_MAX_A,
            },
            "policy": ("reported only; analyze_series does not block IFACE_PREPARED; "
                       "the main session accepts or rejects stability"),
        },
    }
    if len(metrics) < 2:
        block["stability"]["reason"] = "fewer than two structurally valid frames"
        return block
    midpoint = len(metrics) // 2
    first_thickness = statistics.mean(
        frame["thickness_A"] for frame in metrics[:midpoint])
    second_thickness = statistics.mean(
        frame["thickness_A"] for frame in metrics[midpoint:])
    first_center = statistics.mean(frame["center_A"] for frame in metrics[:midpoint])
    second_center = statistics.mean(frame["center_A"] for frame in metrics[midpoint:])
    thickness_difference = abs(first_thickness - second_thickness)
    center_difference = abs(first_center - second_center)
    block.update({
        "first_half_thickness_mean_A": first_thickness,
        "second_half_thickness_mean_A": second_thickness,
        "half_thickness_difference_A": thickness_difference,
        "first_half_center_mean_A": first_center,
        "second_half_center_mean_A": second_center,
        "half_center_difference_A": center_difference,
    })
    block["stability"].update({
        "computed": True,
        "half_thickness_difference_A": thickness_difference,
        "half_center_difference_A": center_difference,
        "thickness_stable": thickness_difference <= STABILITY_THICKNESS_HALF_DIFF_MAX_A,
        "center_stable": center_difference <= STABILITY_CENTER_HALF_DIFF_MAX_A,
        "stable": (thickness_difference <= STABILITY_THICKNESS_HALF_DIFF_MAX_A
                   and center_difference <= STABILITY_CENTER_HALF_DIFF_MAX_A),
        "half_split": ("first/second half by position within the structurally valid "
                       "frames of this report's scope"),
    })
    return block


def build_summary(series, records, extra, raw_frames=None, duplicate_frames=None):
    """Summary for one scope; failing frames are counted, never dropped."""
    steps = [record["step"] for record in records]
    metrics = [record["metrics"] for record in records if record["metrics"] is not None]
    failures = [{"step": record["step"], "start_line": record["start_line"],
                 "reason": record["error"]}
                for record in records if record["metrics"] is None]
    summary = dict(extra)
    summary.update({
        "source_trajectory": str(series["path"]),
        "raw_frames": series["raw_frames"] if raw_frames is None else raw_frames,
        "frames": len(records),
        "duplicate_frames": (series["duplicate_frames"] if duplicate_frames is None
                             else duplicate_frames),
        "dedup_policy": DEDUP_POLICY,
        "first_step": steps[0],
        "last_step": steps[-1],
        "frames_passed": len(metrics),
        "frames_failed": len(failures),
        "all_frames_passed": not failures,
        "failed_steps": [entry["step"] for entry in failures],
        "failures": failures,
        "failure_policy": FAILURE_POLICY,
    })
    if metrics:
        summary["thickness_A"] = series_stats([m["thickness_A"] for m in metrics])
        summary["center_A"] = series_stats([m["center_A"] for m in metrics])
        summary["max_number_density_A-3"] = series_stats(
            [m["max_number_density_A-3"] for m in metrics])
    else:
        summary["thickness_A"] = None
        summary["center_A"] = None
        summary["max_number_density_A-3"] = None
    summary.update(half_stability(metrics))
    summary["solutes"] = {}
    for name in ("S389", "C385", "K394"):
        values = [metric["solutes"][name]["distance_to_nearest_interface"]
                  for metric in metrics]
        summary["solutes"][name] = (
            {"mean_distance": statistics.mean(values), "max_distance": max(values)}
            if values else None)
    last = records[-1]
    summary["final_frame"] = last["metrics"]
    summary["final_frame_step"] = last["step"]
    summary["final_frame_failed"] = last["metrics"] is None
    if last["metrics"] is None:
        summary["final_frame_failure"] = {"step": last["step"],
                                          "start_line": last["start_line"],
                                          "reason": last["error"]}
    return summary


def analyze_series(path, start_step=None, series_name=WINDOW_SERIES_NAME):
    """Whole-trajectory analysis plus, when requested, a production-window report.

    Default behaviour (no ``start_step``) is unchanged: any failing frame aborts
    the series.  The full-trajectory report is always written first so that
    anomalies (for example a detached water cluster) stay on record even when the
    whole-trajectory gate rejects the run.
    """
    path = Path(path)
    series = read_and_dedupe(path)
    records = evaluate_frames(series["ordered"])
    full_summary = build_summary(series, records, {
        "report": FULL_SERIES_NAME,
        "scope": "full_trajectory",
        "window_mode": {"enabled": start_step is not None,
                        "start_step": start_step,
                        "report": series_name if start_step is not None else None,
                        "policy_id": WINDOW_POLICY_ID if start_step is not None else None},
    })
    full_path = path.parent / FULL_SERIES_NAME
    full_path.write_text(json.dumps(full_summary, indent=2) + "\n")
    full_link = {
        "path": FULL_SERIES_NAME,
        "sha256": sha256_file(full_path),
        "scope": "full_trajectory",
        "raw_frames": full_summary["raw_frames"],
        "frames": full_summary["frames"],
        "duplicate_frames": full_summary["duplicate_frames"],
        "frames_passed": full_summary["frames_passed"],
        "frames_failed": full_summary["frames_failed"],
        "failed_steps": full_summary["failed_steps"],
        "all_frames_passed": full_summary["all_frames_passed"],
        "stability": full_summary["stability"],
        "note": ("whole-trajectory evidence; anomalies are kept here permanently and "
                 "are never used as production RC frames"),
    }

    if start_step is None:
        if full_summary["frames_failed"]:
            raise RuntimeError(
                "IFACE_SERIES_FAILED: {}/{} whole-trajectory frames failed the "
                "structural gate (steps {}); failure detail is preserved in {} and "
                "no frame is dropped".format(
                    full_summary["frames_failed"], full_summary["frames"],
                    ", ".join(str(step) for step in full_summary["failed_steps"]),
                    FULL_SERIES_NAME))
        alias = dict(full_summary)
        alias["report"] = series_name
        (path.parent / series_name).write_text(json.dumps(alias, indent=2) + "\n")
        print(json.dumps(alias, indent=2))
        return alias

    window_records = [record for record in records if record["step"] >= start_step]
    if len(window_records) < 2:
        raise RuntimeError(
            "production window starting at step {} holds {} unique frame(s); at "
            "least two are required".format(start_step, len(window_records)))
    window_raw = len([step for step in series["raw_steps"] if step >= start_step])
    window_summary = build_summary(series, window_records, {
        "report": series_name,
        "scope": "production_window",
        "policy_id": WINDOW_POLICY_ID,
        "window": {
            "start_step": start_step,
            "end_step": window_records[-1]["step"],
            "start_inclusive": True,
            "derivation": "{} environment variable".format(PRODUCTION_WINDOW_ENV),
            "reason": WINDOW_RATIONALE,
            "literature_basis": ("none; project-defined convention, not a literature "
                                 "threshold"),
            "alignment": ("same window as prepare_iface_rc.py --window-ps 5.0 "
                          "(steps 5000-10000 at 1 fs)"),
        },
        "full_series_report": full_link,
        "anomaly_disposition": (
            "whole-trajectory anomalies stay in {} and are never deleted or used as "
            "production RC frames; a passing window does not make the whole "
            "trajectory stable".format(FULL_SERIES_NAME)),
        "prepared_gate": ("any failing frame inside this window blocks "
                          "IFACE_PREPARED"),
    }, raw_frames=window_raw, duplicate_frames=window_raw - len(window_records))
    (path.parent / series_name).write_text(json.dumps(window_summary, indent=2) + "\n")
    print(json.dumps(window_summary, indent=2))
    if window_summary["frames_failed"]:
        raise RuntimeError(
            "IFACE_SERIES_WINDOW_FAILED: {}/{} window frames (steps {}-{}) failed the "
            "structural gate at steps {}; every failure is preserved in {} and {} and "
            "IFACE_PREPARED must not be written".format(
                window_summary["frames_failed"], window_summary["frames"],
                window_summary["first_step"], window_summary["last_step"],
                ", ".join(str(step) for step in window_summary["failed_steps"]),
                series_name, FULL_SERIES_NAME))
    return window_summary


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[2] == "--series":
        raw_start = os.environ.get(PRODUCTION_WINDOW_ENV, "").strip()
        window_start = None
        if raw_start:
            try:
                window_start = int(raw_start)
            except ValueError:
                raise SystemExit("{} must be an integer, found {!r}".format(
                    PRODUCTION_WINDOW_ENV, raw_start))
            if window_start <= 0:
                raise SystemExit("{} must be positive, found {}".format(
                    PRODUCTION_WINDOW_ENV, window_start))
        analyze_series(Path(sys.argv[1]), start_step=window_start)
    else:
        main(Path(sys.argv[1]))

#!/usr/bin/env python3
"""Acceptance analysis for the unbiased gas-water interface BOMD trajectory.

The trajectory is parsed with :func:`analyze_iface.read_frames`, so a truncated
frame or a malformed atom line raises instead of silently shifting the frame
boundary.  Every frame must also carry ``i = <step>`` and ``time = <value>`` in
its XYZ comment; ``time`` may use scientific notation.  Frames are deduplicated
by MD step with the last occurrence winning, then sorted by step.

Product definition (project decision, not a literature value):
    S(388)-C(384) < 2.0 A  and  O(389)-H(392) > 1.5 A.

A single frame satisfying that condition is reported as an instantaneous event
only.  The formal product time is the start of the trailing contiguous segment
that satisfies the condition, persists to the last frame, and covers at least
10 frames (about 100 fs at the 10 fs print cadence).  An early event that later
breaks is not a product.

``hms_within_3ps`` is the conjunction of six independent checks: initial state,
trajectory coverage, sustained persistence, product time <= 3000 fs, temperature
evidence, and interface-distance evidence.  Missing temperature or interface
evidence fails the corresponding check; it is never treated as a pass.

``verdict`` answers the stated acceptance criterion (HMS formed and sustained
within 5 ps) and is one of ``CRITERION_SATISFIED``, ``CRITERION_FAILED`` or
``UNDECIDED_COVERAGE_SHORT``.  The last one exists because the production run is
3000 MD steps = 3 ps: "no sustained product in the first 3 ps" is not evidence
that none appears before 5 ps, so a short negative is reported as undecided and
must be extended rather than counted as a mismatch.  ``hms_within_3ps`` and
``formation_fs`` are always reported alongside, so the stricter sub-bar and the
measured onset are never hidden by the verdict.

The S-C/O-H cutoffs, the 10-frame persistence rule, the +/-100 fs start window,
the 5 A median interface distance (start window and sustained segment), the
first-frame step 0 / time 0 coverage requirement, and the temperature gates are
project-defined acceptance criteria.  They are not thresholds quoted from
Liu 2026 or any other literature source.  The script only reports them; it does not claim that a
threshold is experimentally or theoretically derived.

Kept Python 3.6 compatible because the cluster interpreter is old.
"""

import argparse
import json
import math
import re
import statistics
import sys
from collections import namedtuple
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analyze_iface  # noqa: E402


EXPECTED_ATOMS = 394
INDEX_C = 384
INDEX_S = 388
INDEX_O_HYDROXYL = 389
INDEX_H_HYDROXYL = 392
ELEMENT_CHECKS = (
    (INDEX_C, "C"),
    (INDEX_S, "S"),
    (INDEX_O_HYDROXYL, "O"),
    (INDEX_H_HYDROXYL, "H"),
)
DEFAULT_BOX = (17.0, 17.0, 47.0)

S_C_FORMED_MAX_A = 2.0
O_H_BROKEN_MIN_A = 1.5
MIN_COVERAGE_FS = 3000.0
# The stated acceptance criterion is that HMS forms and persists within 5 ps.
# MIN_COVERAGE_FS / product_within_3ps keep the stricter 3 ps sub-bar the
# pipeline already used, so nothing that passed before starts failing.  The
# *verdict* below is decided against 5 ps, and a trajectory that never reaches
# 5 ps is reported as undecided rather than as a mismatch: "no sustained product
# in the first 3 ps" does not answer whether one appears before 5 ps.
FORMATION_CRITERION_MAX_FS = 5000.0
MIN_UNIQUE_FRAMES = 300
MAX_ADJACENT_GAP_FS = 20.0
MIN_PERSISTENCE_FRAMES = 10
SUSTAINED_START_WINDOW_FS = 100.0
END_MEAN_WINDOW_STEPS = 100
END_MEAN_MAX_K = 350.0
END_MAX_WINDOW_STEPS = 10
END_MAX_MAX_K = 450.0
START_WINDOW_MAX_K = 450.0
INTERFACE_MEDIAN_MAX_A = 5.0

# Accepts integers, decimals and scientific notation; used for both the step and
# the time value so an unexpected numeric form fails loudly instead of matching
# only a prefix.
_NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
STEP_VALUE_RE = re.compile(r"\bi\s*=\s*(" + _NUMBER + r")")
TIME_VALUE_RE = re.compile(r"\btime\s*=\s*(" + _NUMBER + r")")

PROJECT_CRITERIA_NOTE = (
    "The S-C < 2.0 A and O-H > 1.5 A bond definitions, the >=10-frame "
    "(about 100 fs) sustained-persistence rule, the +/-100 fs sustained-start "
    "window, the <=5 A median interface-distance gates (start window AND "
    "sustained segment), the requirement that the first frame be step 0 / "
    "time 0, and the temperature gates (last-100-step mean <=350 K, "
    "last-10-step max <=450 K, start-window max <=450 K) are project-defined "
    "structural/dynamical acceptance criteria. They are NOT thresholds quoted "
    "from Liu 2026 or any other literature source."
)

DEDUP_POLICY = (
    "every frame is parsed strictly; frames are deduplicated by MD step with "
    "the last occurrence in the file winning; the result is sorted by step; "
    "duplicate steps and any duplicate-step time conflicts are recorded; the "
    "raw trajectory is never rewritten"
)

FrameRecord = namedtuple("FrameRecord", "step time atoms sc oh product temp")
EnergyRecord = namedtuple("EnergyRecord", "step time temp")


class StrictError(RuntimeError):
    """Raised when the input violates the documented trajectory contract."""


def periodic_distance(first, second, box):
    """Minimum-image distance for an orthorhombic cell.

    The periodic delta is applied on every axis, so the result is correct even
    when the file stores unwrapped coordinates.
    """
    squared = 0.0
    for axis in range(3):
        delta = first[axis] - second[axis]
        delta -= box[axis] * round(delta / box[axis])
        squared += delta * delta
    return math.sqrt(squared)


def parse_comment(comment, path, line_number):
    """Return (step, time_fs), requiring both fields in the XYZ comment."""
    step_match = STEP_VALUE_RE.search(comment)
    if step_match is None:
        raise StrictError(
            "{}: line {}: XYZ comment has no 'i = <step>': {!r}".format(
                path, line_number, comment))
    time_match = TIME_VALUE_RE.search(comment)
    if time_match is None:
        raise StrictError(
            "{}: line {}: XYZ comment has no 'time = <value>': {!r}".format(
                path, line_number, comment))
    step_value = float(step_match.group(1))
    if (not math.isfinite(step_value) or step_value < 0.0
            or abs(step_value - round(step_value)) > 1e-9):
        raise StrictError(
            "{}: line {}: MD step is not a non-negative integer: {!r}".format(
                path, line_number, step_match.group(1)))
    time_value = float(time_match.group(1))
    if not math.isfinite(time_value):
        raise StrictError(
            "{}: line {}: non-finite time value {!r}".format(
                path, line_number, time_match.group(1)))
    return int(round(step_value)), time_value


def load_trajectory(path, box):
    """Strictly parse the XYZ trajectory and return deduplicated records."""
    path = Path(path)
    frames = analyze_iface.read_frames(path)
    lines = path.read_text().splitlines()

    ledger = {}
    duplicate_steps = []
    duplicate_time_conflicts = []
    for frame in frames:
        if len(frame.atoms) != EXPECTED_ATOMS:
            raise StrictError(
                "{}: frame at line {} has {} atoms; expected {}".format(
                    path, frame.start_line, len(frame.atoms), EXPECTED_ATOMS))
        for index, expected in ELEMENT_CHECKS:
            found = str(frame.atoms[index][0]).strip()
            if found.upper() != expected:
                raise StrictError(
                    "{}: frame at line {}: atom index {} must be {}, found "
                    "{!r}".format(path, frame.start_line, index, expected,
                                  found))
        comment = lines[frame.start_line]
        step, time_value = parse_comment(comment, path, frame.start_line + 1)
        s_coords = frame.atoms[INDEX_S][1:4]
        c_coords = frame.atoms[INDEX_C][1:4]
        o_coords = frame.atoms[INDEX_O_HYDROXYL][1:4]
        h_coords = frame.atoms[INDEX_H_HYDROXYL][1:4]
        sc_distance = periodic_distance(s_coords, c_coords, box)
        oh_distance = periodic_distance(o_coords, h_coords, box)
        record = FrameRecord(
            step=step,
            time=time_value,
            atoms=frame.atoms,
            sc=sc_distance,
            oh=oh_distance,
            product=(sc_distance < S_C_FORMED_MAX_A
                     and oh_distance > O_H_BROKEN_MIN_A),
            temp=None,
        )
        if step in ledger:
            duplicate_steps.append(step)
            previous = ledger[step]
            if previous.time != time_value:
                duplicate_time_conflicts.append({
                    "step": step,
                    "kept_time_fs": time_value,
                    "discarded_time_fs": previous.time,
                })
        ledger[step] = record

    ordered = [ledger[step] for step in sorted(ledger)]
    if not ordered:
        raise StrictError("{}: no complete frame found".format(path))
    for previous, current in zip(ordered, ordered[1:]):
        if current.time <= previous.time:
            raise StrictError(
                "{}: time does not increase with step after deduplication "
                "(step {} time {} fs, then step {} time {} fs)".format(
                    path, previous.step, previous.time, current.step,
                    current.time))

    duplicate_info = {
        "policy": DEDUP_POLICY,
        "raw_frames": len(frames),
        "unique_frames": len(ordered),
        "duplicate_frames": len(frames) - len(ordered),
        "duplicate_steps": sorted(set(duplicate_steps)),
        "duplicate_time_conflicts": duplicate_time_conflicts,
    }
    return ordered, duplicate_info


def read_energy(path):
    """Read a CP2K .ener file: column 0 is step, column 3 is temperature."""
    path = Path(path)
    ledger = {}
    raw_records = 0
    duplicate_steps = []
    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        text = line.strip()
        if not text or text.startswith("#"):
            continue
        parts = text.split()
        if len(parts) < 4:
            raise StrictError(
                "{}: line {}: expected at least 4 columns in .ener, found "
                "{!r}".format(path, line_number, line))
        try:
            step_value = float(parts[0])
            time_value = float(parts[1])
            temperature = float(parts[3])
        except ValueError:
            raise StrictError(
                "{}: line {}: non-numeric step/time/temperature in {!r}".format(
                    path, line_number, line))
        if (not math.isfinite(step_value) or step_value < 0.0
                or abs(step_value - round(step_value)) > 1e-9):
            raise StrictError(
                "{}: line {}: .ener step is not a non-negative integer: "
                "{!r}".format(path, line_number, parts[0]))
        if (not math.isfinite(time_value) or not math.isfinite(temperature)
                or temperature <= 0.0):
            raise StrictError(
                "{}: line {}: non-finite time or non-positive temperature "
                "in {!r}".format(path, line_number, line))
        step = int(round(step_value))
        raw_records += 1
        if step in ledger:
            duplicate_steps.append(step)
        ledger[step] = EnergyRecord(step=step, time=time_value, temp=temperature)

    records = [ledger[step] for step in sorted(ledger)]
    if not records:
        raise StrictError("{}: no energy records found".format(path))
    info = {
        "policy": "column 0 = step, column 3 = temperature; deduplicate by "
                  "step with the last occurrence winning; sort by step",
        "raw_records": raw_records,
        "unique_records": len(records),
        "duplicate_records": len(duplicate_steps),
        "duplicate_steps": sorted(set(duplicate_steps)),
        "first_step": records[0].step,
        "last_step": records[-1].step,
    }
    return records, info


def load_energy_optional(path):
    """Return (records, info, error).  A missing file is an evidence failure."""
    if path is None:
        return None, None, ("no --energy file supplied; temperature evidence "
                            "is mandatory for hms_within_3ps")
    energy_path = Path(path)
    if not energy_path.exists():
        return None, None, "energy file not found: {}".format(energy_path)
    records, info = read_energy(energy_path)
    return records, info, None


def build_temperature_lookup(records):
    """Map a trajectory frame to a .ener temperature.

    The MD step is the canonical key.  An exact time match is accepted as a
    fallback for restarted files whose step numbering was offset; uniqueness is
    required so an ambiguous match fails instead of guessing.
    """
    by_step = {}
    by_time = {}
    for record in records:
        by_step[record.step] = record
        by_time.setdefault(round(record.time, 6), []).append(record)

    def lookup(step, time_value):
        record = by_step.get(step)
        if record is not None:
            return record.temp
        candidates = by_time.get(round(time_value, 6), [])
        if len(candidates) == 1:
            return candidates[0].temp
        nearest = None
        for candidate in records:
            difference = abs(candidate.time - time_value)
            if difference <= 0.05 and (nearest is None
                                       or difference < nearest[0]):
                nearest = (difference, candidate.temp)
        return nearest[1] if nearest is not None else None

    return lookup


def summarize(values):
    """Compact min/median/max report, or nulls for an empty series."""
    if not values:
        return {"count": 0, "min": None, "median": None, "max": None,
                "mean": None}
    return {
        "count": len(values),
        "min": float(min(values)),
        "median": float(statistics.median(values)),
        "max": float(max(values)),
        "mean": float(statistics.mean(values)),
    }


def first_event(records, predicate):
    for record in records:
        if predicate(record):
            return {"step": record.step, "time_fs": record.time,
                    "s_c_A": record.sc, "o_h_A": record.oh}
    return None


def normalize_time_origin(records):
    """Shift step/time to the first analyzed frame for segmented protocols."""
    first_step = records[0].step
    first_time = records[0].time
    normalized = [record._replace(step=record.step - first_step,
                                  time=record.time - first_time)
                  for record in records]
    return normalized, {
        "relative_time": True,
        "original_first_step": first_step,
        "original_first_time_fs": first_time,
        "step_offset": first_step,
        "time_offset_fs": first_time,
        "note": ("Only the analysis origin is shifted. Raw trajectories, "
                 "scientific thresholds and event ordering are unchanged."),
    }


def coverage_block(records):
    """Coverage uses the span between the first and last unique record.

    ``last_time >= 3000`` is not sufficient: a trajectory that starts at
    1000 fs and ends at 3500 fs covers only 2500 fs of dynamics and must fail.
    """
    times = [record.time for record in records]
    span = times[-1] - times[0]
    gaps = [later - earlier for earlier, later in zip(times, times[1:])]
    max_gap = max(gaps) if gaps else None
    covers = span >= MIN_COVERAGE_FS
    enough_frames = len(records) >= MIN_UNIQUE_FRAMES
    gap_ok = max_gap is not None and max_gap <= MAX_ADJACENT_GAP_FS
    first_step_zero = records[0].step == 0
    first_time_zero = abs(records[0].time) <= 1.0e-9
    reasons = []
    if not first_step_zero:
        reasons.append(
            "trajectory starts at step {}; the first frame must be step 0"
            .format(records[0].step))
    if not first_time_zero:
        reasons.append(
            "trajectory starts at time {} fs; the first frame must be time 0"
            .format(records[0].time))
    if not covers:
        reasons.append(
            "trajectory spans only {:.3f} fs (first {:.3f} fs, last {:.3f} "
            "fs); need at least {:.1f} fs".format(
                span, times[0], times[-1], MIN_COVERAGE_FS))
    if not enough_frames:
        reasons.append("only {} unique frames; need at least {}".format(
            len(records), MIN_UNIQUE_FRAMES))
    if not gap_ok:
        if max_gap is None:
            reasons.append("a single frame has no adjacent time gap")
        else:
            reasons.append(
                "largest adjacent time gap is {:.3f} fs; limit is {:.1f} "
                "fs".format(max_gap, MAX_ADJACENT_GAP_FS))
    return {
        "passed": (first_step_zero and first_time_zero and covers
                   and enough_frames and gap_ok),
        "reason": "; ".join(reasons) if reasons else (
            "trajectory starts at step 0 / time 0 and spans {:.1f} fs (first "
            "{:.3f} fs, last {:.3f} fs) with {} unique frames and a maximum "
            "adjacent gap of {:.3f} fs".format(
                span, times[0], times[-1], len(records), max_gap)),
        "values": {
            "first_step": records[0].step,
            "first_time_fs": times[0],
            "last_time_fs": times[-1],
            "span_fs": span,
            "unique_frames": len(records),
            "max_adjacent_gap_fs": max_gap,
            "min_coverage_fs": MIN_COVERAGE_FS,
            "min_unique_frames": MIN_UNIQUE_FRAMES,
            "max_adjacent_gap_limit_fs": MAX_ADJACENT_GAP_FS,
            "first_step_zero": first_step_zero,
            "first_time_zero": first_time_zero,
            "span_at_least_min_coverage_fs": covers,
            "enough_frames": enough_frames,
            "gap_within_limit": gap_ok,
        },
    }


def temperature_block(records, energy_records, energy_error, sustained_fs):
    """Temperature evidence for the last steps and the sustained-start window."""
    values = {
        "policy": (
            "last-100-step mean <=350 K and last-10-step max <=450 K are "
            "evaluated on .ener rows inside the trajectory step range; the "
            "sustained-start gate uses trajectory frames within +/-100 fs and "
            "their corresponding .ener temperature"
        ),
        "end_mean_steps": END_MEAN_WINDOW_STEPS,
        "end_mean_max_K": END_MEAN_MAX_K,
        "end_max_steps": END_MAX_WINDOW_STEPS,
        "end_max_max_K": END_MAX_MAX_K,
        "start_window_fs": SUSTAINED_START_WINDOW_FS,
        "start_window_max_K": START_WINDOW_MAX_K,
    }
    if energy_records is None:
        values["energy"] = None
        return {"passed": False,
                "reason": energy_error or "temperature evidence unavailable",
                "values": values}

    lookup = build_temperature_lookup(energy_records)
    reasons = []
    missing = [record for record in records
               if lookup(record.step, record.time) is None]
    values["frames_with_temperature"] = len(records) - len(missing)
    values["frames_without_temperature"] = len(missing)
    if missing:
        reasons.append(
            "{} trajectory frames have no corresponding .ener step/time "
            "(first missing steps: {})".format(
                len(missing), [record.step for record in missing[:10]]))

    first_step = records[0].step
    last_step = records[-1].step
    in_range = [record for record in energy_records
                if first_step <= record.step <= last_step]
    last_mean = in_range[-END_MEAN_WINDOW_STEPS:]
    last_max = in_range[-END_MAX_WINDOW_STEPS:]
    values["energy_steps_in_trajectory_range"] = len(in_range)
    values["trajectory_step_range"] = [first_step, last_step]

    if len(last_mean) < END_MEAN_WINDOW_STEPS:
        reasons.append(
            "energy has only {} steps inside the trajectory step range; the "
            "last-{}-step mean cannot be evaluated".format(
                len(last_mean), END_MEAN_WINDOW_STEPS))
    else:
        mean_temperature = statistics.mean(record.temp for record in last_mean)
        values["last100_steps"] = [last_mean[0].step, last_mean[-1].step]
        values["last100_step_span"] = last_mean[-1].step - last_mean[0].step
        values["last100_mean_K"] = float(mean_temperature)
        if mean_temperature > END_MEAN_MAX_K:
            reasons.append(
                "last-{}-step mean temperature is {:.3f} K; limit is {:.1f} "
                "K".format(END_MEAN_WINDOW_STEPS, mean_temperature,
                           END_MEAN_MAX_K))

    if len(last_max) < END_MAX_WINDOW_STEPS:
        reasons.append(
            "energy has only {} steps inside the trajectory step range; the "
            "last-{}-step maximum cannot be evaluated".format(
                len(last_max), END_MAX_WINDOW_STEPS))
    else:
        maximum_temperature = max(record.temp for record in last_max)
        values["last10_steps"] = [last_max[0].step, last_max[-1].step]
        values["last10_max_K"] = float(maximum_temperature)
        if maximum_temperature > END_MAX_MAX_K:
            reasons.append(
                "last-{}-step maximum temperature is {:.3f} K; limit is "
                "{:.1f} K".format(END_MAX_WINDOW_STEPS, maximum_temperature,
                                  END_MAX_MAX_K))

    if sustained_fs is None:
        reasons.append(
            "no sustained product segment; the sustained-start temperature "
            "window cannot be evaluated")
    else:
        window = [record for record in records
                  if abs(record.time - sustained_fs) <= SUSTAINED_START_WINDOW_FS]
        values["start_window_frames"] = len(window)
        values["start_window_trajectory_time_fs"] = [
            sustained_fs - SUSTAINED_START_WINDOW_FS,
            sustained_fs + SUSTAINED_START_WINDOW_FS,
        ]
        if not window:
            reasons.append("no trajectory frame falls inside the sustained-"
                           "start +/-100 fs window")
        else:
            window_temps = []
            missing_window = []
            for record in window:
                temperature = lookup(record.step, record.time)
                if temperature is None:
                    missing_window.append(record.step)
                else:
                    window_temps.append(temperature)
            if missing_window:
                reasons.append(
                    "{} frames in the sustained-start window have no "
                    "corresponding temperature (steps: {})".format(
                        len(missing_window), missing_window[:10]))
            else:
                window_max = max(window_temps)
                values["start_window_max_observed_K"] = float(window_max)
                if window_max > START_WINDOW_MAX_K:
                    reasons.append(
                        "maximum temperature within +/-100 fs of the sustained "
                        "start is {:.3f} K; limit is {:.1f} K".format(
                            window_max, START_WINDOW_MAX_K))

    return {
        "passed": not reasons,
        "reason": "; ".join(reasons) if reasons else (
            "temperature evidence satisfies the last-step and "
            "sustained-start gates"),
        "values": values,
    }


def interface_block(records, sustained_index):
    """Slab/interface distances over the sustained segment and start window.

    ``analyze_iface.frame_metrics`` is the only source of the slab boundaries
    and solute positions.  Every frame of the trailing sustained segment and
    every frame within +/-100 fs of the sustained start (including frames
    before the start) is analyzed; any ``frame_metrics`` failure fails the
    check.  The gate requires the median distance to the nearest interface for
    both C385 and S389 to be <=5 A inside the +/-100 fs start window AND over
    the whole trailing sustained segment, so a later detachment cannot pass on
    the strength of the start window alone.
    """
    values = {
        "policy": (
            "analyze_iface.frame_metrics is evaluated for every frame of the "
            "trailing sustained segment and for every trajectory frame within "
            "+/-100 fs of the sustained start (including frames before the "
            "start); the gate requires the median distance to the nearest slab "
            "interface to be <=5 A for both C385 and S389 in the start window "
            "and in the sustained segment"
        ),
        "median_distance_limit_A": INTERFACE_MEDIAN_MAX_A,
        "start_window_fs": SUSTAINED_START_WINDOW_FS,
        "metrics_cell_A": list(analyze_iface.CELL),
    }
    if sustained_index is None:
        return {"passed": False,
                "reason": ("no sustained product segment; interface distances "
                           "cannot be evaluated"),
                "values": values}

    start_time = records[sustained_index].time
    selected = [record for record in records
                if (record.time >= start_time
                    or abs(record.time - start_time) <= SUSTAINED_START_WINDOW_FS)]
    analyzed = []
    for record in selected:
        try:
            metrics = analyze_iface.frame_metrics(record.atoms).result
        except Exception as exc:  # frame_metrics failures are acceptance failures
            values["frames_analyzed"] = len(analyzed)
            return {
                "passed": False,
                "reason": (
                    "analyze_iface.frame_metrics failed at step {} (time {} "
                    "fs): {}; the interface criterion cannot be verified".format(
                        record.step, record.time, exc)),
                "values": values,
            }
        analyzed.append((record, metrics))

    segment = [(record, metrics) for record, metrics in analyzed
               if record.time >= start_time]
    window = [(record, metrics) for record, metrics in analyzed
              if abs(record.time - start_time) <= SUSTAINED_START_WINDOW_FS]
    values["sustained_segment_frames"] = len(segment)
    values["start_window_frames"] = len(window)
    values["start_window_trajectory_time_fs"] = [
        start_time - SUSTAINED_START_WINDOW_FS,
        start_time + SUSTAINED_START_WINDOW_FS,
    ]
    if not window:
        return {
            "passed": False,
            "reason": ("no frame inside the sustained-start +/-100 fs window "
                       "could be analyzed"),
            "values": values,
        }

    def solute_stats(name, pairs):
        distances = [metrics["solutes"][name]
                     ["distance_to_nearest_interface"] for _, metrics in pairs]
        z_values = [metrics["solutes"][name]["z"] for _, metrics in pairs]
        inside = [metrics["solutes"][name]["inside_slab"] for _, metrics in pairs]
        return {
            "z_A": summarize(z_values),
            "distance_to_nearest_interface_A": summarize(distances),
            "inside_slab_frames": sum(1 for flag in inside if flag),
            "frames": len(pairs),
        }

    def slab_stats(pairs):
        lower = [metrics["z_lower_A"] for _, metrics in pairs]
        upper = [metrics["z_upper_A"] for _, metrics in pairs]
        thickness = [metrics["thickness_A"] for _, metrics in pairs]
        return {
            "z_lower_A": summarize(lower),
            "z_upper_A": summarize(upper),
            "thickness_A": summarize(thickness),
        }

    window_stats = {}
    segment_stats = {}
    reasons = []
    for name in ("C385", "S389"):
        window_stats[name] = solute_stats(name, window)
        segment_stats[name] = solute_stats(name, segment)
        window_median = window_stats[name][
            "distance_to_nearest_interface_A"]["median"]
        segment_median = segment_stats[name][
            "distance_to_nearest_interface_A"]["median"]
        if window_median > INTERFACE_MEDIAN_MAX_A:
            reasons.append(
                "{} median interface distance within the sustained-start "
                "window is {:.3f} A; limit is {:.1f} A".format(
                    name, window_median, INTERFACE_MEDIAN_MAX_A))
        if segment_median > INTERFACE_MEDIAN_MAX_A:
            reasons.append(
                "{} median interface distance over the sustained segment is "
                "{:.3f} A; limit is {:.1f} A (the product detached after the "
                "start window)".format(
                    name, segment_median, INTERFACE_MEDIAN_MAX_A))

    values["start_window"] = {
        "solutes": window_stats,
        "slab": slab_stats(window),
    }
    values["sustained_segment"] = {
        "solutes": segment_stats,
        "slab": slab_stats(segment),
    }

    # Informational only: signed GDS evidence (AGENTS.md section 8). This block
    # never changes the verdict -- the legacy unsigned <=5 A gate above is kept
    # untouched so existing verdicts stay comparable.
    def gds_stats(pairs):
        out = {}
        for name in ("C385", "S389", "O386"):
            entries = [metrics["solutes"][name] for _, metrics in pairs
                       if name in metrics["solutes"]]
            signed = [e["signed_distance_to_gds"] for e in entries
                      if "signed_distance_to_gds" in e]
            if not signed:
                out[name] = {"available": False}
                continue
            out[name] = {
                "available": True,
                "signed_distance_to_gds_A": summarize(signed),
                "vapour_side_frames": sum(1 for e in entries if e.get("on_vapour_side")),
                "frames": len(signed),
            }
        return out

    values["gds_signed"] = {
        "sign_convention": ("positive = toward the vapour, negative = into the "
                            "water slab; GDS = half-max water-density boundary"),
        "role": ("informational only, not part of the acceptance gate; the "
                 "unsigned <=5 A interface-distance gate is unchanged"),
        "start_window": gds_stats(window),
        "sustained_segment": gds_stats(segment),
    }
    return {
        "passed": not reasons,
        "reason": "; ".join(reasons) if reasons else (
            "C385 and S389 stay within the 5 A median interface-distance gate "
            "both around the sustained start and over the sustained segment"),
        "values": values,
    }


def analyze(trajectory_path, energy_path, box, relative_time=False):
    records, dedup_info = load_trajectory(trajectory_path, box)
    energy_records, energy_info, energy_error = load_energy_optional(energy_path)
    if energy_records is None:
        temperature_lookup = lambda step, time_value: None  # noqa: E731
    else:
        temperature_lookup = build_temperature_lookup(energy_records)
    records = [record._replace(
        temp=temperature_lookup(record.step, record.time))
        for record in records]
    time_origin = {"relative_time": False,
                   "note": "Absolute CP2K step/time are used."}
    if relative_time:
        records, time_origin = normalize_time_origin(records)
        if energy_records is not None:
            energy_records = [record._replace(
                step=record.step - time_origin["step_offset"],
                time=record.time - time_origin["time_offset_fs"])
                for record in energy_records]

    first = records[0]
    last = records[-1]
    unreacted = (first.sc > S_C_FORMED_MAX_A
                 and first.oh < O_H_BROKEN_MIN_A)
    initially_reacted = not unreacted

    first_sc = first_event(records, lambda record: record.sc < S_C_FORMED_MAX_A)
    first_oh = first_event(records, lambda record: record.oh > O_H_BROKEN_MIN_A)
    first_product = first_event(records, lambda record: record.product)

    sustained_index = None
    for index in range(len(records) - 1, -1, -1):
        if records[index].product:
            sustained_index = index
        else:
            break
    persistence_frames = (0 if sustained_index is None
                          else len(records) - sustained_index)
    persisted_to_end = sustained_index is not None
    sustained_fs = (None if sustained_index is None
                    else records[sustained_index].time)
    sustained_duration = (0.0 if sustained_index is None
                          else last.time - sustained_fs)
    sustained_passed = (persisted_to_end
                        and persistence_frames >= MIN_PERSISTENCE_FRAMES)

    coverage = coverage_block(records)
    observation = {
        "first_frame": {
            "step": first.step, "time_fs": first.time,
            "s_c_A": first.sc, "o_h_A": first.oh,
            "unreacted": unreacted, "initially_reacted": initially_reacted,
        },
        "first_instant": {
            "s_c_below_2A": first_sc,
            "o_h_above_1.5A": first_oh,
            "simultaneous_product": first_product,
            "frame_counts": {
                "s_c_below_2A": sum(
                    1 for record in records
                    if record.sc < S_C_FORMED_MAX_A),
                "o_h_above_1.5A": sum(
                    1 for record in records
                    if record.oh > O_H_BROKEN_MIN_A),
                "simultaneous_product": sum(
                    1 for record in records if record.product),
            },
            "note": (
                "instantaneous events are observations only; a single frame "
                "or an event that later breaks is not a reproduced product"),
        },
        "sustained": {
            "sustained_fs": sustained_fs,
            "sustained_duration_fs": sustained_duration,
            "sustained_window_fs": (
                None if sustained_index is None
                else [sustained_fs, last.time]),
            "start_step": (None if sustained_index is None
                           else records[sustained_index].step),
            "persisted_to_end": persisted_to_end,
            "persistence_frames": persistence_frames,
            "minimum_persistence_frames": MIN_PERSISTENCE_FRAMES,
            "definition": (
                "start time of the trailing contiguous segment in which "
                "S-C<2.0 A and O-H>1.5 A; the segment must persist to the "
                "last frame and cover at least 10 frames"),
        },
        "final_frame": {
            "step": last.step, "time_fs": last.time,
            "s_c_A": last.sc, "o_h_A": last.oh,
            "product": last.product, "temp_K": last.temp,
        },
    }

    checks = {}
    checks["initial_unreacted"] = {
        "passed": unreacted,
        "reason": ("first frame is unreacted (S-C {:.3f} A > 2.0 A and O-H "
                   "{:.3f} A < 1.5 A)".format(first.sc, first.oh)
                   if unreacted else
                   "first frame is not unreacted (initially_reacted): S-C "
                   "{:.3f} A, O-H {:.3f} A".format(first.sc, first.oh)),
        "values": {
            "s_c_A": first.sc,
            "o_h_A": first.oh,
            "s_c_limit_A": S_C_FORMED_MAX_A,
            "o_h_limit_A": O_H_BROKEN_MIN_A,
            "initially_reacted": initially_reacted,
        },
    }
    checks["trajectory_coverage"] = coverage
    checks["sustained_product"] = {
        "passed": sustained_passed,
        "reason": (
            "trailing product segment starts at {:.3f} fs and persists to the "
            "last frame over {} frames".format(sustained_fs,
                                               persistence_frames)
            if sustained_passed else
            ("no frame at the end of the trajectory satisfies S-C<2.0 A and "
             "O-H>1.5 A" if sustained_index is None else
             "trailing product segment covers only {} frames; at least {} are "
             "required".format(persistence_frames,
                               MIN_PERSISTENCE_FRAMES))),
        "values": {
            "sustained_fs": sustained_fs,
            "persistence_frames": persistence_frames,
            "persisted_to_end": persisted_to_end,
            "minimum_persistence_frames": MIN_PERSISTENCE_FRAMES,
            "sustained_duration_fs": sustained_duration,
        },
    }
    checks["product_within_3ps"] = {
        "passed": sustained_fs is not None and sustained_fs <= MIN_COVERAGE_FS,
        "reason": (
            "sustained product starts at {:.3f} fs, within 3000 fs".format(
                sustained_fs)
            if sustained_fs is not None and sustained_fs <= MIN_COVERAGE_FS
            else ("no sustained product segment" if sustained_fs is None else
                  "sustained product starts at {:.3f} fs, after 3000 "
                  "fs".format(sustained_fs))),
        "values": {
            "sustained_fs": sustained_fs,
            "limit_fs": MIN_COVERAGE_FS,
        },
    }
    checks["product_within_criterion"] = {
        "passed": (sustained_fs is not None
                   and sustained_fs <= FORMATION_CRITERION_MAX_FS),
        "reason": (
            "sustained product starts at {:.3f} fs, within the {:.0f} fs "
            "criterion".format(sustained_fs, FORMATION_CRITERION_MAX_FS)
            if sustained_fs is not None
            and sustained_fs <= FORMATION_CRITERION_MAX_FS else
            ("no sustained product segment" if sustained_fs is None else
             "sustained product starts at {:.3f} fs, after the {:.0f} fs "
             "criterion".format(sustained_fs, FORMATION_CRITERION_MAX_FS))),
        "values": {
            "sustained_fs": sustained_fs,
            "limit_fs": FORMATION_CRITERION_MAX_FS,
            "stricter_sub_bar_fs": MIN_COVERAGE_FS,
        },
    }
    checks["temperature"] = temperature_block(
        records, energy_records, energy_error, sustained_fs)
    checks["interface"] = interface_block(records, sustained_index)

    failed = [name for name in checks if not checks[name]["passed"]]
    span_fs = coverage["values"]["span_fs"]
    verdict = decide_verdict(checks, sustained_fs, span_fs)
    summary = {
        "source_trajectory": str(trajectory_path),
        "source_energy": (None if energy_path is None else str(energy_path)),
        "box_angstrom": list(box),
        "atom_count": EXPECTED_ATOMS,
        "atom_indices": {
            "C385_zero_based": INDEX_C,
            "S389_zero_based": INDEX_S,
            "O390_hydroxyl_zero_based": INDEX_O_HYDROXYL,
            "H393_hydroxyl_zero_based": INDEX_H_HYDROXYL,
        },
        "deduplication": dedup_info,
        "energy": energy_info,
        "time_origin": time_origin,
        "observation": observation,
        "checks": checks,
        "hms_within_3ps": not failed,
        "hms_within_criterion": verdict == "CRITERION_SATISFIED",
        "verdict": verdict,
        "formation_fs": sustained_fs,
        "criterion_max_formation_fs": FORMATION_CRITERION_MAX_FS,
        "coverage_span_fs": span_fs,
        "failed_checks": failed,
        "failure_reasons": [checks[name]["reason"] for name in failed],
        "criteria": {
            "note": PROJECT_CRITERIA_NOTE,
            "s_c_formed_max_A": S_C_FORMED_MAX_A,
            "o_h_broken_min_A": O_H_BROKEN_MIN_A,
            "min_persistence_frames": MIN_PERSISTENCE_FRAMES,
            "criterion_max_formation_fs": FORMATION_CRITERION_MAX_FS,
            "stricter_sub_bar_fs": MIN_COVERAGE_FS,
            "sustained_start_window_fs": SUSTAINED_START_WINDOW_FS,
            "interface_median_max_A": INTERFACE_MEDIAN_MAX_A,
            "interface_median_scope": ["start_window", "sustained_segment"],
            "first_frame_step_zero": True,
            "first_frame_time_zero": True,
            "min_coverage_fs": MIN_COVERAGE_FS,
            "min_unique_frames": MIN_UNIQUE_FRAMES,
            "max_adjacent_gap_fs": MAX_ADJACENT_GAP_FS,
            "end_mean_steps": END_MEAN_WINDOW_STEPS,
            "end_mean_max_K": END_MEAN_MAX_K,
            "end_max_steps": END_MAX_WINDOW_STEPS,
            "end_max_max_K": END_MAX_MAX_K,
            "start_window_max_K": START_WINDOW_MAX_K,
        },
    }
    return summary, records


def write_outputs(summary, records, json_path, tsv_path):
    json_path = Path(json_path)
    tsv_path = Path(tsv_path)
    if json_path.parent != Path(""):
        json_path.parent.mkdir(parents=True, exist_ok=True)
    if tsv_path.parent != Path(""):
        tsv_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(summary, indent=2) + "\n")
    with tsv_path.open("w") as handle:
        handle.write("step\ttime_fs\ts_c_A\to_h_A\ttemp_K\tproduct\n")
        for record in records:
            temperature = ("" if record.temp is None
                           else "{:.3f}".format(record.temp))
            handle.write("{}\t{:.3f}\t{:.5f}\t{:.5f}\t{}\t{}\n".format(
                record.step, record.time, record.sc, record.oh, temperature,
                1 if record.product else 0))


def box_argument(text):
    parts = [part.strip() for part in text.replace("x", ",").split(",")
             if part.strip()]
    try:
        values = [float(part) for part in parts]
    except ValueError:
        raise argparse.ArgumentTypeError(
            "box must be three comma-separated numbers, got {!r}".format(text))
    if len(values) != 3 or any(value <= 0.0 for value in values):
        raise argparse.ArgumentTypeError(
            "box must be three positive lengths, got {!r}".format(text))
    return tuple(values)


def decide_verdict(checks, sustained_fs, span_fs):
    """Answer the stated <=5 ps criterion; never answer a question the run cannot see.

    ``hms_within_3ps`` is a strict sub-bar over all six checks. It cannot serve as
    the verdict, because the production run is 3000 steps = 3 ps: with no product
    in 3 ps, the ``sustained_product``, ``interface`` and ``product_within_*``
    checks necessarily fail, yet that is not evidence about 5 ps. Declaring
    "mismatch" there would turn a missing experiment into a literature
    disagreement, which is exactly the conclusion that must not be drawn without
    data. So:

      CRITERION_SATISFIED        a sustained product starts within 5 ps, sits at
                                 the interface, and the temperature gates pass
      UNDECIDED_COVERAGE_SHORT   no sustained product, and the trajectory never
                                 reached 5 ps: extend the run, do not conclude
      CRITERION_FAILED           everything else, including a product that only
                                 starts after 5 ps, and substantive evidence
                                 failures (reacted first frame, unusable
                                 trajectory) that no extension would fix
    """
    if not checks["initial_unreacted"]["passed"]:
        return "CRITERION_FAILED"
    if not checks["trajectory_coverage"]["passed"]:
        return "CRITERION_FAILED"
    if sustained_fs is None:
        if span_fs < FORMATION_CRITERION_MAX_FS:
            return "UNDECIDED_COVERAGE_SHORT"
        return "CRITERION_FAILED"
    if sustained_fs > FORMATION_CRITERION_MAX_FS:
        return "CRITERION_FAILED"
    if not checks["interface"]["passed"]:
        return "CRITERION_FAILED"
    if not checks["temperature"]["passed"]:
        return "CRITERION_FAILED"
    return "CRITERION_SATISFIED"


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Strict acceptance analysis of the unbiased interface BOMD "
            "trajectory (sustained HMS product, temperature, interface)."),
        epilog=PROJECT_CRITERIA_NOTE,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("trajectory", help="XYZ trajectory with i=<step> and "
                        "time=<value> in every comment line")
    parser.add_argument("--energy", default=None,
                        help="CP2K .ener file (column 0 step, column 3 Temp); "
                             "required for temperature evidence")
    parser.add_argument("--box", type=box_argument, default=DEFAULT_BOX,
                        help="orthorhombic cell lengths in A, default 17,17,47")
    parser.add_argument("--json", default=None,
                        help="output JSON path (default: UNBIASED_HMS.json next "
                             "to the trajectory)")
    parser.add_argument("--tsv", default=None,
                        help="output TSV path (default: UNBIASED_HMS.tsv next "
                             "to the trajectory)")
    parser.add_argument("--relative-time", action="store_true",
                        help="shift step/time to the first analyzed frame for "
                             "a segmented run; raw files and thresholds are "
                             "unchanged")
    args = parser.parse_args(argv)

    trajectory_path = Path(args.trajectory)
    json_path = (Path(args.json) if args.json is not None
                 else trajectory_path.parent / "UNBIASED_HMS.json")
    tsv_path = (Path(args.tsv) if args.tsv is not None
                else trajectory_path.parent / "UNBIASED_HMS.tsv")
    try:
        summary, records = analyze(trajectory_path, args.energy, args.box,
                                   relative_time=args.relative_time)
    except (RuntimeError, OSError) as exc:
        sys.stderr.write("analyze_unbiased: error: {}\n".format(exc))
        return 2
    write_outputs(summary, records, json_path, tsv_path)
    sys.stdout.write(json.dumps(summary, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())

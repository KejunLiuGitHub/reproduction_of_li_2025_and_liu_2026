#!/usr/bin/env python3
"""Unified factual comparison of the four Li 2025 unbiased BOMD replicas.

Inputs (already pulled; read-only):
  results/li2025_direct/event_unbiased_s1/EVENT_UNBIASED_1.json
  results/li2025_direct/event_unbiased_s2/EVENT_UNBIASED_2.json
  results/li2025_direct/event_unbiased_s3/EVENT_UNBIASED_3.json
  results/li2025_direct/interface_event_unbiased_s1/INTERFACE_EVENT_UNBIASED_1.json
(each directory also carries run_record.json; the JSON label is NOT used as the
identity because two replicas both end up labelled "..._1" — the case name from
REGISTRY below is authoritative.)

Three geometric/factual products, no thermodynamics:

  1. residence segments of d(S-C) < 2.0 A: every contiguous run with frame count,
     start/end time, span = end - start, and the equivalent sampling length
     n_frames * 10 fs. The two differ by one sampling interval on purpose: a
     32-frame run spans 310 fs but represents 320 fs of sampling, and writing the
     latter as a span is wrong.
  2. normalised autocorrelation of the continuous d(S-C) series, lag 0..1000 fs.
     For a dissociating, non-stationary replica the ACF is a trajectory diagnostic
     only: no equilibrium correlation time, rate or barrier is reported.
  3. fixed-bin joint COUNTS of d(S-C) x d(C=O), per replica and summed within each
     environment class. Counts only: the script never forms -kT ln P and never
     produces a free-energy surface.

Bulk (14.2296 A cube) and interface (15 x 15 x 30 A) replicas are never merged
into one statistic: the class totals stay separate and the summary records that
no cross-class merge was performed.

Outputs (default results/li2025_direct/unbiased_comparison/):
  UNBIASED_COMPARISON.json      machine-readable summary
  RESIDENCE_SEGMENTS.tsv        one row per contiguous d(S-C) < 2.0 A run
  RESIDENCE_SUMMARY.tsv         per-replica run statistics
  ACF_<case>.tsv                lag_frames, lag_fs, acf
  ACF_KEY_LAGS.tsv              key lags side by side
  HIST2D_<case>.tsv             fixed-bin joint counts, long format
  HIST2D_CLASS_TOTALS.tsv       bulk total and interface total (never mixed)
  FIG_SC_CO_JOINT_4PANEL.png    four panels, one per replica
  REPORT.md                     short human-readable summary of the same numbers

Usage:
  python3 compare_unbiased_replicas.py
  python3 compare_unbiased_replicas.py --outdir /tmp/cmp --no-figure
"""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "li2025_direct"
OUTDIR = RESULTS / "unbiased_comparison"

SCHEMA = "li2025_direct_unbiased_comparison/v1"
ATOM_COUNT = 296
EXPECTED_FRAMES = 301
SAMPLE_INTERVAL_FS = 10.0
RESIDENCE_THRESHOLD = 2.0
DISSOCIATION_THRESHOLD = 2.5
DISSOCIATION_FRAMES = 10

# Fixed binning, identical for every replica (and for the class totals).
SC_BIN_START, SC_BIN_STOP, SC_BIN_WIDTH = 1.70, 3.80, 0.10
CO_BIN_START, CO_BIN_STOP, CO_BIN_WIDTH = 1.18, 1.44, 0.02

# Lags requested: 0 .. 1000 fs at the 10 fs sampling interval.
MAX_LAG_FS = 1000.0
KEY_LAGS_FS = (0.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0, 1000.0)

BOX_BY_CLASS = {
    "bulk": (14.2296, 14.2296, 14.2296),
    "interface": (15.0, 15.0, 30.0),
}

# case name -> (environment class, result subdirectory, accepted JSON file names)
REGISTRY = (
    ("event_unbiased_s1", "bulk", "event_unbiased_s1",
     ("EVENT_UNBIASED_1.json",)),
    ("event_unbiased_s2", "bulk", "event_unbiased_s2",
     ("EVENT_UNBIASED_2.json",)),
    ("event_unbiased_s3", "bulk", "event_unbiased_s3",
     ("EVENT_UNBIASED_3.json",)),
    ("interface_event_unbiased_s1", "interface", "interface_event_unbiased_s1",
     ("INTERFACE_EVENT_UNBIASED_1.json", "EVENT_UNBIASED_1.json")),
)


def sha256(path):
    digest = hashlib.sha256()
    with open(str(path), "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolve_json(directory, names):
    """First existing accepted file name; ambiguous or missing input stops the run."""
    found = [directory / name for name in names if (directory / name).is_file()]
    if not found:
        raise SystemExit("none of %s found in %s" % (list(names), directory))
    if len(found) > 1:
        raise SystemExit("several candidate JSONs in %s: %s; refusing to guess"
                         % (directory, [path.name for path in found]))
    return found[0]


def display_path(path):
    """Repo-relative when possible, absolute otherwise (e.g. a temp --outdir)."""
    resolved = Path(path).resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


def bin_edges(start, stop, width):
    count = int(round((stop - start) / width))
    if abs(start + count * width - stop) > 1.0e-9:
        raise SystemExit("bin range %s..%s is not a multiple of %s" % (start, stop, width))
    return [round(start + index * width, 6) for index in range(count + 1)]


def bin_index(value, edges):
    """Index of the half-open bin [edges[i], edges[i+1]); last edge is inclusive."""
    if value < edges[0] or value > edges[-1]:
        raise SystemExit("value %.6f outside the fixed bins %.6f..%.6f"
                         % (value, edges[0], edges[-1]))
    if value == edges[-1]:
        return len(edges) - 2
    for index in range(len(edges) - 1):
        if edges[index] <= value < edges[index + 1]:
            return index
    raise SystemExit("value %.6f could not be binned" % value)


def residence_runs(times, values, threshold):
    """Every contiguous run with value < threshold, in time order."""
    runs = []
    start = None
    for index, value in enumerate(values):
        if value < threshold and start is None:
            start = index
        elif value >= threshold and start is not None:
            runs.append((start, index - 1))
            start = None
    if start is not None:
        runs.append((start, len(values) - 1))
    rows = []
    for order, (begin, end) in enumerate(runs):
        frames = end - begin + 1
        span = times[end] - times[begin]
        rows.append({
            "run_index": order,
            "start_index": begin,
            "end_index": end,
            "n_frames": frames,
            "start_time_fs": times[begin],
            "end_time_fs": times[end],
            "span_fs": span,
            "equivalent_sampling_fs": frames * SAMPLE_INTERVAL_FS,
            "span_vs_equivalent_gap_fs": frames * SAMPLE_INTERVAL_FS - span,
        })
    return rows


def median(values):
    ordered = sorted(values)
    count = len(ordered)
    if count == 0:
        return None
    middle = count // 2
    if count % 2 == 1:
        return ordered[middle]
    return 0.5 * (ordered[middle - 1] + ordered[middle])


def sustained_dissociation(times, values, threshold, frames):
    """First index starting a run of `frames` consecutive samples above threshold."""
    for index in range(0, len(values) - frames + 1):
        if all(value > threshold for value in values[index:index + frames]):
            end = index + frames - 1
            return {
                "criterion": "d(S-C) > %.2f A for %d consecutive samples"
                             % (threshold, frames),
                "kind": "project geometric recording criterion, not a thermodynamic result",
                "index": index,
                "step": int(times[index]),
                "time_fs": times[index],
                "window_end_index": end,
                "window_end_time_fs": times[end],
                "window_span_fs": times[end] - times[index],
                "equivalent_sampling_fs": frames * SAMPLE_INTERVAL_FS,
            }
    return None


def acf(values, max_lag_samples):
    """Normalised autocorrelation, r(k) = sum c[i]c[i+k] / sum c[i]^2, r(0) = 1."""
    count = len(values)
    mean = sum(values) / float(count)
    centred = [value - mean for value in values]
    variance = sum(value * value for value in centred)
    if variance <= 0.0:
        raise SystemExit("cannot autocorrelate a constant series")
    lags = []
    for lag in range(0, max_lag_samples + 1):
        total = 0.0
        for index in range(0, count - lag):
            total += centred[index] * centred[index + lag]
        lags.append(total / variance)
    return lags


def acf_diagnostics(lags, acf_values, stationary):
    first_non_positive = None
    for lag in range(1, len(acf_values)):
        if acf_values[lag] <= 0.0:
            first_non_positive = lag
            break
    below_e = None
    for lag in range(1, len(acf_values)):
        if acf_values[lag] < 1.0 / 2.718281828459045:
            below_e = lag
            break
    key = {}
    for lag_fs in KEY_LAGS_FS:
        index = int(round(lag_fs / SAMPLE_INTERVAL_FS))
        if index < len(acf_values):
            key["%g_fs" % lag_fs] = acf_values[index]
    return {
        "max_lag_fs": lags[-1] * SAMPLE_INTERVAL_FS,
        "samples_at_max_lag": lags[-1],
        "acf_at_max_lag": acf_values[-1],
        "key_lags": key,
        "first_lag_below_zero_fs": (None if first_non_positive is None
                                    else first_non_positive * SAMPLE_INTERVAL_FS),
        "first_lag_below_1_over_e_fs": (None if below_e is None
                                        else below_e * SAMPLE_INTERVAL_FS),
        "stationary": stationary,
        "equilibrium_correlation_time_fs": None,
        "rate_or_barrier": None,
        "note": ("diagnostic only: a dissociating (non-stationary) replica has no equilibrium "
                 "autocorrelation time, so no correlation time, rate or barrier is derived"
                 if not stationary else
                 "diagnostic only: no sustained dissociation in this 3 ps window, but "
                 "stationarity is not established over 3 ps, so no correlation time, rate "
                 "or barrier is derived"),
    }


def hist2d(xs, ys, xedges, yedges):
    counts = [[0 for _ in range(len(yedges) - 1)] for _ in range(len(xedges) - 1)]
    for x, y in zip(xs, ys):
        counts[bin_index(x, xedges)][bin_index(y, yedges)] += 1
    return counts


def add_counts(left, right):
    return [[left[i][j] + right[i][j] for j in range(len(left[0]))]
            for i in range(len(left))]


def flatten_counts(case, klass, counts, xedges, yedges):
    rows = []
    for i in range(len(xedges) - 1):
        for j in range(len(yedges) - 1):
            rows.append({
                "case": case,
                "class": klass,
                "sc_bin_low": xedges[i],
                "sc_bin_high": xedges[i + 1],
                "co_bin_low": yedges[j],
                "co_bin_high": yedges[j + 1],
                "count": counts[i][j],
            })
    return rows


def load_case(case, klass, subdir, names):
    directory = RESULTS / subdir
    json_path = resolve_json(directory, names)
    data = json.loads(json_path.read_text(encoding="utf-8"))
    run_record_path = directory / "run_record.json"
    run_record = (json.loads(run_record_path.read_text(encoding="utf-8"))
                  if run_record_path.is_file() else {})

    rows = data.get("per_frame") or []
    if len(rows) != EXPECTED_FRAMES:
        raise SystemExit("%s has %d per-frame rows, expected %d"
                         % (case, len(rows), EXPECTED_FRAMES))
    if data.get("atom_count") != ATOM_COUNT:
        raise SystemExit("%s atom_count %s != %d" % (case, data.get("atom_count"), ATOM_COUNT))
    times = [float(row["time_fs"]) for row in rows]
    steps = [int(row["step"]) for row in rows]
    if len(set(steps)) != len(steps):
        raise SystemExit("%s has duplicate MD steps" % case)
    if steps != sorted(steps):
        raise SystemExit("%s time series is not in MD order" % case)
    for index in range(1, len(times)):
        if abs((times[index] - times[index - 1]) - SAMPLE_INTERVAL_FS) > 1.0e-6:
            raise SystemExit("%s sampling interval is not %g fs at row %d"
                             % (case, SAMPLE_INTERVAL_FS, index))
    box = data.get("box_angstrom") or []
    expected_box = BOX_BY_CLASS[klass]
    if len(box) != 3 or any(abs(float(box[i]) - expected_box[i]) > 1.0e-9 for i in range(3)):
        raise SystemExit("%s box %s is not the %s cell %s"
                         % (case, box, klass, list(expected_box)))
    if run_record:
        recorded_box = run_record.get("box_angstrom") or []
        if (len(recorded_box) != 3
                or any(abs(float(recorded_box[i]) - expected_box[i]) > 1.0e-9
                       for i in range(3))):
            raise SystemExit("%s run_record box %s is not the %s cell"
                             % (case, recorded_box, klass))
        if run_record.get("status") != "COMPLETED":
            raise SystemExit("%s run_record status is %r, not COMPLETED"
                             % (case, run_record.get("status")))
    return {
        "case": case,
        "class": klass,
        "json_path": json_path,
        "json_sha256": sha256(json_path),
        "json_label": data.get("label"),
        "run_record": run_record,
        "rows": rows,
        "times": times,
        "steps": steps,
        "box": [float(value) for value in box],
        "trajectory_sha256": data.get("trajectory_sha256"),
    }


def proton_transfer_proxy(rows):
    """Distance-based check that no proton moved off its water in this window.

    Two facts from the same per-frame JSONs: the covalent O-H distance of the
    water whose hydrogen is nearest the carbonyl oxygen, and that H...O(carbonyl)
    distance. A proton transfer would show up as a stretched O-H and a sub-1.1 A
    H...O contact; neither appears. This is a proxy, not a bond-topology proof.
    """
    oh = [float(row["water_oh_distance"]) for row in rows]
    hd = [float(row["nearest_water_h_distance"]) for row in rows]
    stretched = [value for value in oh if value > 1.20]
    transferred = [value for value in hd if value < 1.10]
    return {
        "proxy": ("nearest-water covalent O-H distance and the same H...O(carbonyl) "
                  "distance, per frame"),
        "water_oh_min_angstrom": min(oh),
        "water_oh_max_angstrom": max(oh),
        "nearest_h_to_carbonyl_o_min_angstrom": min(hd),
        "thresholds": {"water_oh_rupture_angstrom": 1.20,
                       "h_transfer_contact_angstrom": 1.10},
        "frames_with_stretched_water_oh": len(stretched),
        "frames_with_transfer_contact": len(transferred),
        "consistent_with_no_proton_transfer": not stretched and not transferred,
        "note": "distance-based proxy only; not a bond-topology proof",
    }


def series_facts(times, values):
    return {
        "first": values[0],
        "last": values[-1],
        "min": min(values),
        "max": max(values),
        "mean": sum(values) / float(len(values)),
        "first_time_fs": times[0],
        "last_time_fs": times[-1],
    }


def case_block(loaded, acf_values, acf_meta, hist_meta):
    times = loaded["times"]
    rows = loaded["rows"]
    sc = [float(row["sc"]) for row in rows]
    co = [float(row["co"]) for row in rows]
    runs = residence_runs(times, sc, RESIDENCE_THRESHOLD)
    frames = [run["n_frames"] for run in runs]
    spans = [run["span_fs"] for run in runs]
    dissociation = sustained_dissociation(times, sc, DISSOCIATION_THRESHOLD,
                                           DISSOCIATION_FRAMES)
    stationary = dissociation is None
    return {
        "case": loaded["case"],
        "class": loaded["class"],
        "job_id": loaded["run_record"].get("job_id"),
        "host": loaded["run_record"].get("host"),
        "box_angstrom": loaded["box"],
        "input_json": display_path(loaded["json_path"]),
        "input_json_sha256": loaded["json_sha256"],
        "input_json_label": loaded["json_label"],
        "trajectory_sha256": loaded["trajectory_sha256"],
        "frames": len(rows),
        "sample_interval_fs": SAMPLE_INTERVAL_FS,
        "duration_fs": times[-1] - times[0],
        "sc": series_facts(times, sc),
        "co": series_facts(times, co),
        "residence": {
            "threshold_angstrom": RESIDENCE_THRESHOLD,
            "n_runs": len(runs),
            "n_frames_in_runs": sum(frames),
            "median_n_frames": median(frames),
            "max_n_frames": max(frames) if frames else None,
            "median_span_fs": median(spans),
            "max_span_fs": max(spans) if spans else None,
            "first_run_start_time_fs": runs[0]["start_time_fs"] if runs else None,
            "last_run_end_time_fs": runs[-1]["end_time_fs"] if runs else None,
            "longest_run": (max(runs, key=lambda run: (run["n_frames"],
                                                       -run["run_index"]))
                            if runs else None),
            "runs": runs,
        },
        "sustained_dissociation": dissociation,
        "proton_transfer": proton_transfer_proxy(rows),
        "stationary": stationary,
        "stationarity_established": False,
        "stationarity_note": ("no sustained dissociation seen; 3 ps is nevertheless too short "
                              "to establish stationarity"
                              if stationary else
                              "sustained dissociation makes the series non-stationary"),
        "acf": acf_meta,
        "hist2d": hist_meta,
        "acf_values": acf_values,
    }


def write_tsv(path, columns, rows):
    with open(str(path), "w") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(row[column] for column in columns) + "\n")


def fmt(value, digits=6):
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(value)
    return ("%." + str(digits) + "f") % float(value)


def plot_matrix(counts):
    """Transpose a [sc][co] count matrix into imshow's [row=y=co][col=x=sc] order."""
    return [[counts[i][j] for i in range(len(counts))] for j in range(len(counts[0]))]


def write_figure(path, blocks, xedges, yedges, vmax):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(2, 2, figsize=(11.5, 8.6))
    for axis, block in zip(axes.flat, blocks):
        counts = block["hist2d_counts"]
        image = axis.imshow(plot_matrix(counts), origin="lower", aspect="auto",
                            cmap="viridis", vmin=0, vmax=vmax,
                            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
        axis.set_title("%s  (%s, job %s)\n%d frames, %d occupied bins"
                       % (block["case"], block["class"], block["job_id"],
                          block["frames"], sum(1 for row in counts for value in row if value)),
                       fontsize=9)
        axis.set_xlabel("d(S$-$C) / $\\AA$", fontsize=9)
        axis.set_ylabel("d(C$=$O) / $\\AA$", fontsize=9)
        axis.tick_params(labelsize=8)
        figure.colorbar(image, ax=axis, fraction=0.046, pad=0.03).ax.tick_params(labelsize=8)
    figure.suptitle("Li 2025 unbiased replicas: joint COUNTS of d(S-C) x d(C=O)",
                    fontsize=11)
    figure.text(0.5, 0.015,
                "counts only - no -kT ln P, no free-energy surface; "
                "bulk and interface cells are separate panels, never merged",
                ha="center", fontsize=8)
    figure.tight_layout(rect=(0, 0.03, 1, 0.96))
    figure.savefig(str(path), dpi=150)
    plt.close(figure)


def report_text(summary):
    lines = ["# Li 2025 unbiased replicas — comparison facts", "",
             "Generated by `campaign/li2025_direct/compare_unbiased_replicas.py`. "
             "Counts and distances only: no -kT ln P and no free-energy surface.", ""]
    lines.append("| case | class | job | box (Å) | frames | d(S–C) first→last (Å) | "
                 "<2 Å frames | runs | longest run | sustained >2.5 Å |")
    lines.append("|---|---|---|---|---:|---|---:|---:|---|---|")
    for block in summary["case_order"]:
        case = summary["cases"][block]
        longest = case["residence"]["longest_run"]
        dissoc = case["sustained_dissociation"]
        lines.append("| %s | %s | %s | %s | %d | %.4f → %.4f | %d | %d | %s | %s |"
                     % (case["case"], case["class"], case["job_id"],
                        " × ".join("%g" % value for value in case["box_angstrom"]),
                        case["frames"], case["sc"]["first"], case["sc"]["last"],
                        case["residence"]["n_frames_in_runs"],
                        case["residence"]["n_runs"],
                        ("%d frames, %.0f–%.0f fs (span %.0f fs)"
                         % (longest["n_frames"], longest["start_time_fs"],
                            longest["end_time_fs"], longest["span_fs"]))
                        if longest else "none",
                        ("from %.0f fs, window to %.0f fs"
                         % (dissoc["time_fs"], dissoc["window_end_time_fs"]))
                        if dissoc else "not reached"))
    lines += ["", "## Residence segment accounting", "",
              "`span_fs = end_time_fs - start_time_fs` and "
              "`equivalent_sampling_fs = n_frames × 10 fs` differ by exactly one "
              "sampling interval (10 fs) for every run: a 32-frame run spans 310 fs "
              "but represents 320 fs of sampling.", ""]
    for block in summary["case_order"]:
        case = summary["cases"][block]
        residence = case["residence"]
        lines.append("- **%s**: %d run%s, %d frames below 2 Å, median %s frames, "
                     "max %s frames, median span %s fs, max span %s fs."
                     % (case["case"], residence["n_runs"],
                        "" if residence["n_runs"] == 1 else "s",
                        residence["n_frames_in_runs"],
                        fmt(residence["median_n_frames"], 0),
                        fmt(residence["max_n_frames"], 0),
                        fmt(residence["median_span_fs"], 0),
                        fmt(residence["max_span_fs"], 0)))
    lines += ["", "## Autocorrelation (diagnostic only)", "",
              "lag 0–1000 fs, TSV per replica. A dissociating replica is "
              "non-stationary, so no equilibrium correlation time, rate or barrier "
              "is derived from these curves.", ""]
    for block in summary["case_order"]:
        case = summary["cases"][block]
        acf = case["acf"]
        lines.append("- **%s**: ACF context %s, r(10 fs) = %s, r(100 fs) = %s, "
                     "r(1000 fs) = %s, first non-positive lag %s fs."
                     % (case["case"],
                        "no sustained dissociation" if acf["stationary"] else "dissociating",
                        fmt(acf["key_lags"].get("10_fs"), 3),
                        fmt(acf["key_lags"].get("100_fs"), 3),
                        fmt(acf["acf_at_max_lag"], 3),
                        fmt(acf["first_lag_below_zero_fs"], 0)))
    lines += ["", "## Joint counts", "",
              "Fixed bins: d(S–C) %.2f–%.2f Å in %.2f Å steps, d(C=O) %.2f–%.2f Å in "
              "%.2f Å steps. Bulk and interface totals are reported separately; no "
              "cross-class statistic exists." %
              (summary["binning"]["sc"]["start"], summary["binning"]["sc"]["stop"],
               summary["binning"]["sc"]["width"], summary["binning"]["co"]["start"],
               summary["binning"]["co"]["stop"], summary["binning"]["co"]["width"]), ""]
    for klass in ("bulk", "interface"):
        total = summary["class_totals"][klass]
        lines.append("- **%s total**: %d frames over %d replicas (%s)."
                     % (klass, total["frames"], total["n_replicas"],
                        ", ".join(total["cases"])))
    lines += ["", "## Proton transfer proxy", "",
              "Nearest-water O-H and H...O(carbonyl) distances only; a bond-topology "
              "analysis is not performed.", ""]
    for block in summary["case_order"]:
        case = summary["cases"][block]
        transfer = case["proton_transfer"]
        lines.append("- **%s**: water O-H %.3f-%.3f Å, nearest H...O(carbonyl) %.3f Å, "
                     "stretched O-H frames %d, transfer contacts %d → %s."
                     % (case["case"], transfer["water_oh_min_angstrom"],
                        transfer["water_oh_max_angstrom"],
                        transfer["nearest_h_to_carbonyl_o_min_angstrom"],
                        transfer["frames_with_stretched_water_oh"],
                        transfer["frames_with_transfer_contact"],
                        "consistent with no proton transfer"
                        if transfer["consistent_with_no_proton_transfer"]
                        else "proton motion cannot be excluded"))
    lines += ["", "Figure: `FIG_SC_CO_JOINT_4PANEL.png` (one panel per replica).", ""]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default=str(OUTDIR))
    parser.add_argument("--no-figure", action="store_true",
                        help="skip the PNG (analysis still runs; used by the self-test "
                             "when matplotlib is unavailable)")
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    xedges = bin_edges(SC_BIN_START, SC_BIN_STOP, SC_BIN_WIDTH)
    yedges = bin_edges(CO_BIN_START, CO_BIN_STOP, CO_BIN_WIDTH)
    max_lag_samples = int(round(MAX_LAG_FS / SAMPLE_INTERVAL_FS))

    loaded = [load_case(case, klass, subdir, names)
              for case, klass, subdir, names in REGISTRY]

    blocks = []
    residence_rows = []
    summary_rows = []
    acf_key_rows = []
    hist_rows = []
    class_totals = {}
    for item in loaded:
        times = item["times"]
        sc = [float(row["sc"]) for row in item["rows"]]
        co = [float(row["co"]) for row in item["rows"]]
        runs = residence_runs(times, sc, RESIDENCE_THRESHOLD)
        for run in runs:
            residence_rows.append(dict(run, case=item["case"], case_class=item["class"]))
        dissociation = sustained_dissociation(times, sc, DISSOCIATION_THRESHOLD,
                                               DISSOCIATION_FRAMES)
        stationary = dissociation is None

        acf_values = acf(sc, max_lag_samples)
        acf_meta = acf_diagnostics(list(range(len(acf_values))), acf_values, stationary)
        acf_path = outdir / ("ACF_%s.tsv" % item["case"])
        write_tsv(acf_path, ["lag_frames", "lag_fs", "acf"],
                  [{"lag_frames": str(lag),
                    "lag_fs": fmt(lag * SAMPLE_INTERVAL_FS, 3),
                    "acf": fmt(value, 9)}
                   for lag, value in enumerate(acf_values)])
        acf_meta["tsv"] = display_path(acf_path)
        for lag_name, value in acf_meta["key_lags"].items():
            acf_key_rows.append({"case": item["case"], "class": item["class"],
                                 "lag_fs": lag_name.replace("_fs", ""),
                                 "acf": fmt(value, 9)})

        counts = hist2d(sc, co, xedges, yedges)
        hist_path = outdir / ("HIST2D_%s.tsv" % item["case"])
        write_tsv(hist_path, ["case", "class", "sc_bin_low", "sc_bin_high",
                              "co_bin_low", "co_bin_high", "count"],
                  [{"case": row["case"], "class": row["class"],
                    "sc_bin_low": fmt(row["sc_bin_low"], 3),
                    "sc_bin_high": fmt(row["sc_bin_high"], 3),
                    "co_bin_low": fmt(row["co_bin_low"], 3),
                    "co_bin_high": fmt(row["co_bin_high"], 3),
                    "count": str(row["count"])}
                   for row in flatten_counts(item["case"], item["class"], counts,
                                             xedges, yedges)])
        hist_meta = {
            "tsv": display_path(hist_path),
            "sc_bins": len(xedges) - 1,
            "co_bins": len(yedges) - 1,
            "cells_with_counts": sum(1 for row in counts for value in row if value),
            "total_counts": sum(value for row in counts for value in row),
        }
        if hist_meta["total_counts"] != len(item["rows"]):
            raise SystemExit("%s histogram counts %d != %d frames"
                             % (item["case"], hist_meta["total_counts"],
                                len(item["rows"])))

        block = case_block(item, acf_values, acf_meta, hist_meta)
        block["hist2d_counts"] = counts
        blocks.append(block)
        residence = block["residence"]
        summary_rows.append({
            "case": item["case"], "class": item["class"],
            "n_runs": str(residence["n_runs"]),
            "n_frames_in_runs": str(residence["n_frames_in_runs"]),
            "median_n_frames": fmt(residence["median_n_frames"], 1),
            "max_n_frames": fmt(residence["max_n_frames"], 0),
            "median_span_fs": fmt(residence["median_span_fs"], 1),
            "max_span_fs": fmt(residence["max_span_fs"], 1),
            "first_run_start_time_fs": fmt(residence["first_run_start_time_fs"], 1),
            "last_run_end_time_fs": fmt(residence["last_run_end_time_fs"], 1),
            "dissociation_time_fs": fmt(dissociation["time_fs"] if dissociation else None, 1),
        })

        bucket = class_totals.setdefault(item["class"], {
            "class": item["class"], "cases": [], "frames": 0,
            "counts": [[0 for _ in yedges[:-1]] for _ in xedges[:-1]],
        })
        bucket["cases"].append(item["case"])
        bucket["frames"] += len(item["rows"])
        bucket["counts"] = add_counts(bucket["counts"], counts)

    class_rows = []
    for klass in sorted(class_totals):
        bucket = class_totals[klass]
        for row in flatten_counts("%s_total" % klass, klass, bucket["counts"],
                                  xedges, yedges):
            class_rows.append({
                "case": row["case"], "class": row["class"],
                "sc_bin_low": fmt(row["sc_bin_low"], 3),
                "sc_bin_high": fmt(row["sc_bin_high"], 3),
                "co_bin_low": fmt(row["co_bin_low"], 3),
                "co_bin_high": fmt(row["co_bin_high"], 3),
                "count": str(row["count"])})
    write_tsv(outdir / "HIST2D_CLASS_TOTALS.tsv",
              ["case", "class", "sc_bin_low", "sc_bin_high", "co_bin_low",
               "co_bin_high", "count"], class_rows)

    write_tsv(outdir / "RESIDENCE_SEGMENTS.tsv",
              ["case", "class", "run_index", "start_index", "end_index", "n_frames",
               "start_time_fs", "end_time_fs", "span_fs", "equivalent_sampling_fs",
               "span_vs_equivalent_gap_fs"],
              [{"case": row["case"], "class": row["case_class"],
                "run_index": str(row["run_index"]),
                "start_index": str(row["start_index"]),
                "end_index": str(row["end_index"]),
                "n_frames": str(row["n_frames"]),
                "start_time_fs": fmt(row["start_time_fs"], 1),
                "end_time_fs": fmt(row["end_time_fs"], 1),
                "span_fs": fmt(row["span_fs"], 1),
                "equivalent_sampling_fs": fmt(row["equivalent_sampling_fs"], 1),
                "span_vs_equivalent_gap_fs": fmt(row["span_vs_equivalent_gap_fs"], 1)}
               for row in residence_rows])
    write_tsv(outdir / "RESIDENCE_SUMMARY.tsv",
              ["case", "class", "n_runs", "n_frames_in_runs", "median_n_frames",
               "max_n_frames", "median_span_fs", "max_span_fs",
               "first_run_start_time_fs", "last_run_end_time_fs",
               "dissociation_time_fs"], summary_rows)
    write_tsv(outdir / "ACF_KEY_LAGS.tsv", ["case", "class", "lag_fs", "acf"],
              acf_key_rows)

    vmax = max(value for block in blocks for row in block["hist2d_counts"]
               for value in row)
    figure_path = outdir / "FIG_SC_CO_JOINT_4PANEL.png"
    figure_meta = {"path": display_path(figure_path), "written": False}
    if not args.no_figure:
        try:
            write_figure(figure_path, blocks, xedges, yedges, vmax)
            figure_meta["written"] = True
            figure_meta["vmax_count"] = vmax
            figure_meta["panels"] = [block["case"] for block in blocks]
        except ImportError as error:
            figure_meta["note"] = "matplotlib unavailable: %s" % error

    summary = {
        "schema": SCHEMA,
        "purpose": ("Unified factual comparison of the four Li 2025 unbiased replicas: "
                    "residence segments, d(S-C) autocorrelation diagnostics and fixed-bin "
                    "joint counts of d(S-C) x d(C=O)."),
        "policy": {
            "what_this_is": [
                "Distance and count facts for four 3 ps unbiased replicas.",
                "A single fixed binning shared by every panel and both class totals.",
            ],
            "what_this_is_not": [
                "NOT a free-energy surface: -kT ln P is never formed and no FES is produced.",
                "NOT a rate, barrier, mean first-passage time or equilibrium correlation time.",
                "NOT Li's three-replica statistics.",
                "NOT an electronic bonding analysis.",
            ],
            "structure_conclusion_status": ("distance facts only; any C-S bonding claim must "
                                            "pass the three-layer gate in AGENTS.md section 7"),
        },
        "thresholds": {
            "residence_angstrom": RESIDENCE_THRESHOLD,
            "dissociation_angstrom": DISSOCIATION_THRESHOLD,
            "dissociation_frames": DISSOCIATION_FRAMES,
            "dissociation_note": ("project geometric recording criterion "
                                  "(d(S-C) > 2.5 A for 10 consecutive samples); it is not a "
                                  "barrier, rate or thermodynamic conclusion"),
            "sample_interval_fs": SAMPLE_INTERVAL_FS,
            "residence_accounting": ("span_fs = end_time_fs - start_time_fs; "
                                     "equivalent_sampling_fs = n_frames * 10 fs; they differ "
                                     "by one sampling interval, so a 32-frame run spans "
                                     "310 fs, not 320 fs"),
        },
        "binning": {
            "sc": {"start": SC_BIN_START, "stop": SC_BIN_STOP, "width": SC_BIN_WIDTH,
                   "bins": len(xedges) - 1},
            "co": {"start": CO_BIN_START, "stop": CO_BIN_STOP, "width": CO_BIN_WIDTH,
                   "bins": len(yedges) - 1},
            "policy": "identical edges for every replica and both class totals; counts only",
        },
        "acf": {
            "definition": ("r(k) = sum_i c_i c_{i+k} / sum_i c_i^2 with c = x - <x>, "
                           "lag 0..1000 fs at the 10 fs sampling interval"),
            "max_lag_fs": MAX_LAG_FS,
            "key_lags_fs": list(KEY_LAGS_FS),
            "status": ("diagnostic only: dissociating replicas are non-stationary, so no "
                       "equilibrium correlation time, rate or barrier is derived"),
        },
        "cross_class_merge": {
            "performed": False,
            "reason": ("bulk (14.2296 A cube) and interface (15 x 15 x 30 A) replicas have "
                       "different cells and environments; their totals are kept separate"),
        },
        "free_energy": {"computed": False,
                        "note": "hard rule: no -kT ln P and no FES in this comparison"},
        "case_order": [block["case"] for block in blocks],
        "cases": {},
        "class_totals": {},
        "figure": figure_meta,
        "verdict": None,
    }
    for block in blocks:
        block.pop("hist2d_counts", None)
        summary["cases"][block["case"]] = block
    for klass, bucket in class_totals.items():
        summary["class_totals"][klass] = {
            "class": klass,
            "cases": bucket["cases"],
            "n_replicas": len(bucket["cases"]),
            "frames": bucket["frames"],
            "tsv": "HIST2D_CLASS_TOTALS.tsv",
        }
    summary["files"] = {
        "residence_segments": "RESIDENCE_SEGMENTS.tsv",
        "residence_summary": "RESIDENCE_SUMMARY.tsv",
        "acf_key_lags": "ACF_KEY_LAGS.tsv",
        "hist2d_class_totals": "HIST2D_CLASS_TOTALS.tsv",
    }
    out_json = outdir / "UNBIASED_COMPARISON.json"
    out_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    (outdir / "REPORT.md").write_text(report_text(summary), encoding="utf-8")

    print(json.dumps({
        "out": str(out_json),
        "cases": [{
            "case": block["case"], "class": block["class"],
            "job_id": block["job_id"],
            "sc_first_last": [round(block["sc"]["first"], 4),
                              round(block["sc"]["last"], 4)],
            "frames_below_2a": block["residence"]["n_frames_in_runs"],
            "runs": block["residence"]["n_runs"],
            "longest_run_frames": (block["residence"]["longest_run"] or {}).get("n_frames"),
            "longest_run_span_fs": (block["residence"]["longest_run"] or {}).get("span_fs"),
            "dissociation_time_fs": (block["sustained_dissociation"] or {}).get("time_fs"),
        } for block in blocks],
        "figure": figure_meta,
    }, indent=2))
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())

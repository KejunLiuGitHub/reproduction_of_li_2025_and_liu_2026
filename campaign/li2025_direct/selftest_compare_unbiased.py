#!/usr/bin/env python3
"""Self-test for the four-replica Li 2025 unbiased comparison.

Covers the pure analysis functions with hand-computable synthetic series, then
runs the real comparison end to end on the pulled result JSONs and checks the
files, the recorded facts and the "no free energy" boundary.

Covered:
  * residence runs: contiguous segmentation, frame count, start/end time,
    span = end - start, equivalent sampling = n_frames * 10 fs and the 10 fs gap
    between them (a 32-frame run spans 310 fs, it does not span 320 fs),
  * residence summary: run count, median and max frame count, first/last run ends,
  * sustained dissociation: 10 consecutive samples above 2.5 A, and the negative
    case where a 9-sample excursion must not trigger,
  * autocorrelation: exact hand-computed values for [0,0,1,1], r(0) = 1, constant
    series rejected, key-lag mapping and the non-stationary disclaimer
    (no equilibrium correlation time, no rate, no barrier),
  * fixed-bin histograms: identical edges, counts summing to the frame count,
    last-edge inclusion, out-of-range rejection and the [sc][co] -> [row=y][col=x]
    transpose used by the figure,
  * end-to-end: the four replicas with their case names, job ids, boxes,
    residence/dissociation facts, 101 ACF rows per replica, class totals that keep
    bulk and interface separate, and a written PNG,
  * determinism: two runs produce byte-identical JSON,
  * Python 3.6 syntax/API and `bash -n` for every script in the directory.

Run: python3 selftest_compare_unbiased.py
Exit 0 and a final SELFTEST_COMPARE_UNBIASED_OK line on success.
"""
import sys

sys.dont_write_bytecode = True

import ast
import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))

import compare_unbiased_replicas as cmp

CHECKS = {"count": 0}
EXPECTED_JOBS = {
    "event_unbiased_s1": "8566509",
    "event_unbiased_s2": "8566870",
    "event_unbiased_s3": "8566871",
    "interface_event_unbiased_s1": "8566899",
}


def ok(label, condition, detail=""):
    CHECKS["count"] += 1
    if not condition:
        raise AssertionError("%s%s" % (label, (": " + detail) if detail else ""))
    print("ok %2d - %s" % (CHECKS["count"], label))


def read_text(path):
    return Path(path).read_text(encoding="utf-8")


def run_process(args):
    process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = process.communicate()[0].decode("utf-8", "replace")
    return process.returncode, output


def check_residence_functions():
    times = [index * 10.0 for index in range(12)]
    values = [2.10, 1.90, 1.80, 1.70, 2.50, 2.05, 2.20, 1.50, 1.60, 2.30, 2.05, 2.01]
    runs = cmp.residence_runs(times, values, 2.0)
    ok("two contiguous runs below 2 A are found", len(runs) == 2, str(len(runs)))
    ok("the first run is frames 1-3 (indices 1..3)",
       runs[0]["start_index"] == 1 and runs[0]["end_index"] == 3
       and runs[0]["n_frames"] == 3)
    ok("the first run spans 20 fs but samples 30 fs",
       runs[0]["start_time_fs"] == 10.0 and runs[0]["end_time_fs"] == 30.0
       and runs[0]["span_fs"] == 20.0 and runs[0]["equivalent_sampling_fs"] == 30.0)
    ok("the isolated frame 7-8 is a 2-frame run", runs[1]["n_frames"] == 2
       and runs[1]["start_index"] == 7 and runs[1]["end_index"] == 8)
    ok("every run records span_vs_equivalent_gap_fs = 10",
       all(abs(run["span_vs_equivalent_gap_fs"] - 10.0) < 1e-9 for run in runs))
    ok("a single-frame run has span 0 fs and equivalent 10 fs",
       all(run["equivalent_sampling_fs"] - run["span_fs"] == 10.0 for run in runs))

    # The 32-frame trap: 1410..1720 fs is 310 fs of span, 320 fs of sampling.
    trap_times = [index * 10.0 for index in range(301)]
    trap_values = [2.5] * 301
    for index in range(141, 173):
        trap_values[index] = 1.9
    trap_runs = cmp.residence_runs(trap_times, trap_values, 2.0)
    ok("a 32-frame run reports n_frames = 32", len(trap_runs) == 1
       and trap_runs[0]["n_frames"] == 32)
    ok("that run spans 1410-1720 fs, i.e. 310 fs (not 320)",
       trap_runs[0]["start_time_fs"] == 1410.0 and trap_runs[0]["end_time_fs"] == 1720.0
       and trap_runs[0]["span_fs"] == 310.0)
    ok("its equivalent sampling is 320 fs, 10 fs longer than the span",
       trap_runs[0]["equivalent_sampling_fs"] == 320.0
       and trap_runs[0]["span_vs_equivalent_gap_fs"] == 10.0)

    ok("median of an odd list", cmp.median([3.0, 1.0, 2.0]) == 2.0)
    ok("median of an even list is the midpoint", cmp.median([1.0, 2.0, 3.0, 6.0]) == 2.5)
    ok("median of an empty list is None", cmp.median([]) is None)


def check_dissociation_function():
    times = [index * 10.0 for index in range(40)]
    values = [1.8] * 9 + [2.7] * 9 + [1.8] + [2.7] * 30
    found = cmp.sustained_dissociation(times, values, 2.5, 10)
    ok("a 9-sample excursion does not trigger the 10-sample criterion",
       found is not None and found["index"] == 19,
       "index=%s" % (found and found["index"]))
    ok("the trigger records the first frame of the 10-sample window",
       found["time_fs"] == 190.0 and found["window_end_time_fs"] == 280.0)
    ok("the criterion reports span 90 fs and equivalent sampling 100 fs",
       found["window_span_fs"] == 90.0 and found["equivalent_sampling_fs"] == 100.0)
    ok("the criterion is labelled geometric, not thermodynamic",
       "geometric" in found["kind"] and "not a thermodynamic" in found["kind"])
    ok("a series that never exceeds the threshold has no trigger",
       cmp.sustained_dissociation(times, [2.0] * 40, 2.5, 10) is None)


def check_acf_functions():
    values = [0.0, 0.0, 1.0, 1.0]
    result = cmp.acf(values, 3)
    ok("acf(0) is exactly 1", result[0] == 1.0)
    ok("acf of [0,0,1,1] is 1, 0.25, -0.5, -0.25",
       abs(result[1] - 0.25) < 1e-12 and abs(result[2] + 0.5) < 1e-12
       and abs(result[3] + 0.25) < 1e-12, str(result))

    series = [1.0, 2.0, 1.5, 3.0, 2.2, 1.1]
    auto = cmp.acf(series, len(series) - 1)
    ok("acf(0) is 1 for an arbitrary series", abs(auto[0] - 1.0) < 1e-12)
    ok("acf returns max_lag + 1 values", len(auto) == len(series))
    try:
        cmp.acf([1.0] * 10, 3)
        raise AssertionError("accepted a constant series")
    except SystemExit:
        ok("a constant series is rejected instead of dividing by zero", True)

    lags = list(range(0, 101))
    stationary = cmp.acf_diagnostics(lags, [1.0] + [0.5] * 100, False)
    ok("non-stationary diagnostics state stationarity was not assumed",
       stationary["stationary"] is False)
    ok("non-stationary diagnostics give no correlation time and no rate",
       stationary["equilibrium_correlation_time_fs"] is None
       and stationary["rate_or_barrier"] is None)
    ok("non-stationary diagnostics say diagnostic only",
       "diagnostic only" in stationary["note"]
       and "non-stationary" in stationary["note"])
    ok("the key lag 10 fs is the value at lag index 1",
       abs(stationary["key_lags"]["10_fs"] - 0.5) < 1e-12)
    ok("max lag is 1000 fs over 100 samples",
       stationary["max_lag_fs"] == 1000.0 and stationary["samples_at_max_lag"] == 100)
    ok("first non-positive lag is reported in fs",
       stationary["first_lag_below_zero_fs"] is None)
    crossing = cmp.acf_diagnostics(lags, [1.0] + [0.4] * 5 + [-0.1] * 95, False)
    ok("a zero crossing is reported at the first non-positive lag",
       crossing["first_lag_below_zero_fs"] == 60.0)


def check_histogram_functions():
    xedges = cmp.bin_edges(1.70, 3.80, 0.10)
    yedges = cmp.bin_edges(1.18, 1.44, 0.02)
    ok("the S-C binning has 22 edges / 21 bins ending at 3.80",
       len(xedges) == 22 and len(yedges) == 14
       and xedges[0] == 1.70 and xedges[-1] == 3.80)
    ok("the C=O binning ends at 1.44", yedges[-1] == 1.44)

    xs = [1.75, 1.75, 3.79, 3.80]
    ys = [1.19, 1.19, 1.43, 1.44]
    counts = cmp.hist2d(xs, ys, xedges, yedges)
    ok("the two identical low points share one bin",
       counts[cmp.bin_index(1.75, xedges)][cmp.bin_index(1.19, yedges)] == 2)
    ok("3.80 A lands in the last bin, not outside it",
       cmp.bin_index(3.80, xedges) == len(xedges) - 2)
    ok("counts sum to the number of samples",
       sum(value for row in counts for value in row) == len(xs))
    ok("an out-of-range value is rejected",
       _raises(SystemExit, cmp.bin_index, 3.95, xedges))

    matrix = cmp.plot_matrix(counts)
    ok("the figure matrix is the [row=y=co][col=x=sc] transpose",
       all(matrix[j][i] == counts[i][j]
           for i in range(len(counts)) for j in range(len(counts[0]))))

    doubled = cmp.add_counts(counts, counts)
    ok("class totals add element-wise",
       all(doubled[i][j] == 2 * counts[i][j]
           for i in range(len(counts)) for j in range(len(counts[0]))))


def _raises(exception, function, *args):
    try:
        function(*args)
        return False
    except exception:
        return True


def without_output_paths(cases):
    """Per-replica numbers, with the output-directory-dependent TSV links removed."""
    stripped = {}
    for case, block in cases.items():
        copy = json.loads(json.dumps(block))
        copy["acf"].pop("tsv", None)
        copy["hist2d"].pop("tsv", None)
        stripped[case] = copy
    return stripped


def check_proton_proxy():
    rows = [{"water_oh_distance": 0.99, "nearest_water_h_distance": 1.45},
            {"water_oh_distance": 1.05, "nearest_water_h_distance": 1.38}]
    clean = cmp.proton_transfer_proxy(rows)
    ok("a clean window is flagged consistent with no proton transfer",
       clean["consistent_with_no_proton_transfer"] is True
       and clean["frames_with_stretched_water_oh"] == 0
       and clean["frames_with_transfer_contact"] == 0)
    ok("the proxy records the O-H and H...O ranges",
       abs(clean["water_oh_min_angstrom"] - 0.99) < 1e-12
       and abs(clean["water_oh_max_angstrom"] - 1.05) < 1e-12
       and abs(clean["nearest_h_to_carbonyl_o_min_angstrom"] - 1.38) < 1e-12)
    stretched = cmp.proton_transfer_proxy(
        rows + [{"water_oh_distance": 1.35, "nearest_water_h_distance": 1.40}])
    ok("a stretched water O-H breaks the no-transfer flag",
       stretched["consistent_with_no_proton_transfer"] is False
       and stretched["frames_with_stretched_water_oh"] == 1)
    contacted = cmp.proton_transfer_proxy(
        rows + [{"water_oh_distance": 1.02, "nearest_water_h_distance": 1.05}])
    ok("a sub-1.1 A H...O(carbonyl) contact breaks the no-transfer flag",
       contacted["consistent_with_no_proton_transfer"] is False
       and contacted["frames_with_transfer_contact"] == 1)
    ok("the proxy is labelled a proxy, not a bond-topology proof",
       "proxy" in clean["note"] and "not a bond-topology proof" in clean["note"])


def check_end_to_end():
    scratch = tempfile.mkdtemp(prefix="li2025_compare_selftest_")
    try:
        outdir = Path(scratch) / "cmp"
        code, output = run_process([sys.executable,
                                    str(HERE / "compare_unbiased_replicas.py"),
                                    "--outdir", str(outdir)])
        ok("the comparison script exits 0", code == 0, output.strip()[-400:])
        expected_files = ["UNBIASED_COMPARISON.json", "RESIDENCE_SEGMENTS.tsv",
                          "RESIDENCE_SUMMARY.tsv", "ACF_KEY_LAGS.tsv",
                          "HIST2D_CLASS_TOTALS.tsv", "REPORT.md",
                          "FIG_SC_CO_JOINT_4PANEL.png"]
        missing = [name for name in expected_files if not (outdir / name).is_file()]
        ok("all summary files and the figure exist", not missing, str(missing))

        summary = json.loads(read_text(outdir / "UNBIASED_COMPARISON.json"))
        ok("schema", summary["schema"] == "li2025_direct_unbiased_comparison/v1")
        ok("four replicas in registry order",
           summary["case_order"] == ["event_unbiased_s1", "event_unbiased_s2",
                                     "event_unbiased_s3",
                                     "interface_event_unbiased_s1"])
        ok("every replica carries its case name, job id and box",
           all(summary["cases"][case]["job_id"] == EXPECTED_JOBS[case]
               for case in summary["case_order"]))
        ok("bulk replicas are the 14.2296 A cube",
           all(summary["cases"][case]["box_angstrom"] == [14.2296, 14.2296, 14.2296]
               for case in ("event_unbiased_s1", "event_unbiased_s2",
                            "event_unbiased_s3")))
        ok("the interface replica is the 15 x 15 x 30 A cell",
           summary["cases"]["interface_event_unbiased_s1"]["box_angstrom"]
           == [15.0, 15.0, 30.0])
        ok("each replica has 301 frames over 3000 fs",
           all(summary["cases"][case]["frames"] == 301
               and summary["cases"][case]["duration_fs"] == 3000.0
               for case in summary["case_order"]))

        s1 = summary["cases"]["event_unbiased_s1"]
        ok("bulk s1 keeps the contact: no sustained dissociation",
           s1["sustained_dissociation"] is None and s1["stationary"] is True)
        ok("no replica claims an established stationary ensemble",
           all(summary["cases"][case]["stationarity_established"] is False
               for case in summary["case_order"]))
        ok("bulk s1 has 214 frames below 2 A in 18 runs",
           s1["residence"]["n_frames_in_runs"] == 214
           and s1["residence"]["n_runs"] == 18)
        ok("bulk s1 longest run is 32 frames spanning 310 fs (320 fs sampled)",
           s1["residence"]["longest_run"]["n_frames"] == 32
           and s1["residence"]["longest_run"]["span_fs"] == 310.0
           and s1["residence"]["longest_run"]["equivalent_sampling_fs"] == 320.0)
        ok("bulk s1 longest run starts at 1410 fs and ends at 1720 fs",
           s1["residence"]["longest_run"]["start_time_fs"] == 1410.0
           and s1["residence"]["longest_run"]["end_time_fs"] == 1720.0)

        for case, when in (("event_unbiased_s2", 160.0),
                           ("event_unbiased_s3", 80.0),
                           ("interface_event_unbiased_s1", 180.0)):
            block = summary["cases"][case]
            ok("%s dissociates (10 samples > 2.5 A) at %.0f fs" % (case, when),
               block["sustained_dissociation"] is not None
               and block["sustained_dissociation"]["time_fs"] == when,
               str(block["sustained_dissociation"]))
            ok("%s is therefore flagged non-stationary in the ACF block" % case,
               block["stationary"] is False
               and block["acf"]["equilibrium_correlation_time_fs"] is None)

        ok("class totals keep bulk and interface apart",
           summary["class_totals"]["bulk"]["frames"] == 903
           and summary["class_totals"]["interface"]["frames"] == 301
           and summary["cross_class_merge"]["performed"] is False)
        ok("every replica keeps its own identity (four distinct case names)",
           len(set(summary["case_order"])) == 4
           and summary["cases"]["interface_event_unbiased_s1"]["input_json"]
           .endswith("interface_event_unbiased_s1/INTERFACE_EVENT_UNBIASED_1.json")
           and summary["cases"]["event_unbiased_s1"]["input_json"]
           .endswith("event_unbiased_s1/EVENT_UNBIASED_1.json"))
        for case in summary["case_order"]:
            transfer = summary["cases"][case]["proton_transfer"]
            ok("%s shows no proton-transfer signature" % case,
               transfer["consistent_with_no_proton_transfer"] is True
               and transfer["frames_with_stretched_water_oh"] == 0
               and transfer["frames_with_transfer_contact"] == 0,
               str(transfer["water_oh_max_angstrom"]))
        ok("no free energy is computed",
           summary["free_energy"]["computed"] is False
           and "no -kT ln P" in summary["free_energy"]["note"])
        ok("no verdict is asserted", summary["verdict"] is None)
        ok("thresholds are recorded",
           summary["thresholds"]["residence_angstrom"] == 2.0
           and summary["thresholds"]["dissociation_angstrom"] == 2.5
           and summary["thresholds"]["dissociation_frames"] == 10
           and summary["thresholds"]["sample_interval_fs"] == 10.0)
        ok("the figure lists four panels and was written",
           summary["figure"]["written"] is True
           and len(summary["figure"]["panels"]) == 4)

        for case in summary["case_order"]:
            block = summary["cases"][case]
            acf_rows = read_text(outdir / block["acf"]["tsv"]).splitlines()
            ok("%s ACF TSV has a header plus 101 lags (0..1000 fs)" % case,
               len(acf_rows) == 102 and acf_rows[1].split("\t")[0] == "0"
               and abs(float(acf_rows[1].split("\t")[2]) - 1.0) < 1e-9)
            ok("%s ACF TSV ends at 1000 fs" % case,
               abs(float(acf_rows[-1].split("\t")[1]) - 1000.0) < 1e-9)
            hist_rows = read_text(outdir / block["hist2d"]["tsv"]).splitlines()
            ok("%s histogram TSV has header plus 21 x 13 cells" % case,
               len(hist_rows) == 1 + 21 * 13)
            total = sum(int(row.split("\t")[6]) for row in hist_rows[1:])
            ok("%s joint counts sum to 301 frames" % case, total == 301, str(total))

        residence_rows = read_text(outdir / "RESIDENCE_SEGMENTS.tsv").splitlines()
        ok("the residence TSV has a header and 23 runs (18 + 1 + 3 + 1)",
           len(residence_rows) == 1 + 23, str(len(residence_rows)))
        ok("the residence TSV names both span and equivalent sampling",
           "span_fs" in residence_rows[0]
           and "equivalent_sampling_fs" in residence_rows[0]
           and "span_vs_equivalent_gap_fs" in residence_rows[0])
        gaps = set(row.split("\t")[10] for row in residence_rows[1:])
        ok("every run has a 10 fs gap between span and equivalent sampling",
           gaps == {"10.0"}, str(sorted(gaps)))

        totals_rows = read_text(outdir / "HIST2D_CLASS_TOTALS.tsv").splitlines()
        bulk_total = sum(int(row.split("\t")[6]) for row in totals_rows[1:]
                         if row.split("\t")[1] == "bulk")
        interface_total = sum(int(row.split("\t")[6]) for row in totals_rows[1:]
                              if row.split("\t")[1] == "interface")
        ok("class-total TSV keeps bulk (903) and interface (301) separate",
           bulk_total == 903 and interface_total == 301)

        report = read_text(outdir / "REPORT.md")
        ok("the report states the span/equivalence distinction",
           "310 fs" in report and "-kT ln P" in report)
        ok("the report records the proton-transfer proxy",
           "proton" in report.lower() and "consistent with no proton transfer" in report)
        ok("the figure file is a non-trivial PNG",
           (outdir / "FIG_SC_CO_JOINT_4PANEL.png").stat().st_size > 20000)

        # Determinism: an identical run must reproduce the JSON byte for byte.
        first = read_text(outdir / "UNBIASED_COMPARISON.json")
        code, _ = run_process([sys.executable,
                              str(HERE / "compare_unbiased_replicas.py"),
                              "--outdir", str(outdir)])
        ok("an identical rerun exits 0", code == 0)
        ok("an identical rerun reproduces the summary JSON byte for byte",
           read_text(outdir / "UNBIASED_COMPARISON.json") == first)

        # --no-figure must skip only the PNG and change no number.
        nofig = Path(scratch) / "nofig"
        code, _ = run_process([sys.executable,
                              str(HERE / "compare_unbiased_replicas.py"),
                              "--outdir", str(nofig), "--no-figure"])
        ok("--no-figure exits 0 and writes no PNG",
           code == 0 and not (nofig / "FIG_SC_CO_JOINT_4PANEL.png").exists())
        rerun = json.loads(read_text(nofig / "UNBIASED_COMPARISON.json"))
        ok("--no-figure records that the PNG was not written",
           rerun["figure"]["written"] is False)
        ok("--no-figure leaves every per-replica number unchanged",
           without_output_paths(rerun["cases"])
           == without_output_paths(json.loads(first)["cases"]))
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def check_guards():
    scratch = tempfile.mkdtemp(prefix="li2025_compare_guards_")
    try:
        directory = Path(scratch)
        for name in ("A.json", "B.json"):
            (directory / name).write_text("{}", encoding="utf-8")
        ok("two candidate JSONs are refused instead of picked",
           _raises(SystemExit, cmp.resolve_json, directory, ("A.json", "B.json")))
        ok("a missing input JSON stops the run",
           _raises(SystemExit, cmp.resolve_json, directory, ("C.json",)))
        ok("the accepted interface file name is registered",
           "interface_event_unbiased_s1"
           == [case for case, klass, _subdir, _names in cmp.REGISTRY
               if klass == "interface"][0])
        ok("a box that is not a multiple of the bin width is refused",
           _raises(SystemExit, cmp.bin_edges, 1.70, 3.85, 0.10))
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def check_python36_and_shell():
    py_files = sorted(HERE.glob("*.py"))
    ok("the comparison script and its self-test are present",
       {"compare_unbiased_replicas.py",
        "selftest_compare_unbiased.py"}.issubset(set(path.name for path in py_files)))
    for path in py_files:
        source = read_text(path)
        try:
            ast.parse(source, filename=str(path), feature_version=(3, 6))
            compile(source, str(path), "exec")
        except SyntaxError as error:
            raise AssertionError("%s is not Python 3.6 syntax: %s" % (path.name, error))
        capture_keyword = "capture" + "_output"
        text_keyword = "te" + "xt"
        for token in (capture_keyword, "missing" + "_ok",
                      "remove" + "prefix", "remove" + "suffix"):
            ok("%s has no %s" % (path.name, token), token not in source)
        pattern = (r"subprocess\.(run|Popen|call|check_call|check_output)"
                   r"\([^)]*?\b(" + text_keyword + "|" + capture_keyword + r")\s*=")
        ok("%s has no subprocess text/capture keywords" % path.name,
           re.search(pattern, source, re.DOTALL) is None)
        ok("%s has no walrus operator" % path.name, (":" + "=") not in source)
    for path in sorted(HERE.glob("*.slurm")):
        ok("%s passes bash -n" % path.name,
           subprocess.call(["bash", "-n", str(path)]) == 0)


def check_no_thermodynamics():
    source = read_text(HERE / "compare_unbiased_replicas.py")
    for token in ("boltzmann", "exp(-", "math.exp", "k_b", "kb_t"):
        ok("the comparison script has no %r conversion" % token,
           token not in source.lower())
    normalized = " ".join(source.split())
    ok("the script states the no-FES rule in its docstring",
       "-kT ln P" in normalized and "never produces a free-energy surface" in normalized)
    readme = read_text(HERE / "README.md").lower().replace("−", "-")
    for token in ("compare_unbiased_replicas", "unbiased_comparison",
                  "residence", "autocorrelation", "-kt ln p",
                  "span_vs_equivalent_gap_fs"):
        ok("README documents %r" % token, token in readme)


def main():
    check_residence_functions()
    check_dissociation_function()
    check_acf_functions()
    check_histogram_functions()
    check_guards()
    check_proton_proxy()
    check_end_to_end()
    check_python36_and_shell()
    check_no_thermodynamics()
    print("SELFTEST_COMPARE_UNBIASED_OK: %d checks passed (residence, dissociation, "
          "acf, histograms, guards, end-to-end, py36, boundaries)" % CHECKS["count"])
    return 0


if __name__ == "__main__":
    sys.exit(main())

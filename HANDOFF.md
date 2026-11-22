# Collaboration handoff

**Status:** 22 November 2026

This repository contains our attempts to reproduce selected interfacial-chemistry results from Li et al. (2025) and Liu et al. (2026). We have not reproduced the reported reaction behaviour or obtained converged free-energy barriers. Our results come from reconstructed setups and do not show that either paper is incorrect.

Automatic follow-up was stopped on 20 November. No new calculations were started for this handoff; any continuation needs a separate decision. The downloadable Releases are still dated 19 November. Since then, both Li original-coordinate controls have finished, but their complete files have not yet been added to a new Release. `metadata/HANDOFF_STATUS.json` records their status and hashes.

## 1. Documents and data

Start with the two progress reports, then consult the instructions for the study you are working on:

- `docs/PROGRESS.md` and `docs/LIU2026_PROGRESS.md`: results, coverage and open questions.
- `docs/METHODS.md`: Li atom mapping, model, software and observable definitions.
- `docs/REPRODUCE.md` and `docs/LIU2026_REPRODUCE.md`: replay commands and numerical checks.
- `docs/DATA.md`: Li Release layout and file selection.
- `docs/ADDITIONAL_RESULTS.md`: collection IDs for files outside the public downloads. `REUSE.md` covers attribution and reuse.

The two studies use different atom maps. Li has 296 atoms, SO₃²⁻ + HCHO, 96 waters and BLYP-D3(BJ); Liu has 394 atoms, HOSO₂⁻ + HCHO + K⁺, 128 waters and PBE-D3. Keep their coordinates, bias histories and free-energy analyses separate.

The public package contains the analysis code, case manifests, selected inputs, replay tests and two Releases: `v2026.11.19` for Li and `liu2026-v2026.11.19` for Liu. Some original author inputs, restart files and starting coordinates are unavailable. Where necessary, we reconstructed the setup and recorded the differences.

The current maintainer should also supply:

- Complete files for Li jobs **8572577 and 8572578**, or access to collection `li2025_hbcontrol_20261118/` in the working archive.
- Any additional collection needed for a chosen task, with its manifest and hashes.
- Working-archive or HPC access if new production calculations are planned. Credentials and private notes are transferred separately, never through this repository.
- A fresh scheduler and remote-write check immediately before any submission, cancellation or replacement. The queue status recorded here is historical.

## 2. Check the public package

Use Python 3.12. The original replay tests used Python 3.12.10 and NumPy 2.2.6.

```bash
git clone https://github.com/KejunLiuGitHub/reproduction_of_li_2025_and_liu_2026.git
cd reproduction_of_li_2025_and_liu_2026
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt

python tools/fetch_data.py --dataset li2025
python tools/fetch_data.py --dataset liu2026
python tools/fetch_data.py --dataset all --verify-only
```

Run one replay from each study and then the portable tests:

```bash
python tools/reproduce_case.py \
  --case natural_original40 --out analysis/li_natural_original40
python tools/reproduce_liu_case.py \
  --case frame_ic_1 --out analysis/liu_frame_ic_1

python tools/test_portable_tools.py
python tools/test_liu_tools.py --case frame_ic_1 -v
```

Use a new output directory for each replay. Previous test results are in `examples/PORTABLE_TEST_REPORT.json`, `examples/LIU2026_REPLAY_TEST_REPORT.json` and `examples/DOWNLOAD_TEST_REPORT.json`.

These commands recalculate descriptors from saved trajectories; they do not rerun CP2K or reconstruct a free-energy surface. The portable tests also run without CP2K. `docs/REPRODUCE.md` explains how to prepare a fresh Li CP2K input.

## 3. Where the results stand

| Study | Observed | Still open |
|---|---|---|
| Li 2025 | Short C–S distances under metadynamics; two selected configurations retained predominantly short distances for 5 ps after bias removal | Formation from separated reactants without bias, a converged barrier and interfacial acceleration |
| Liu 2026 | Short distances under guidance; retention after controlled proton placement in two prepared arms | Formation from the intended reactant complex without bias, the reported 0.6 kcal/mol barrier and converged sampling |

A trajectory that starts with a short C–S distance tests retention, not formation from separated reactants. Guided event times are not unbiased reaction times. The five Liu `frame_ic_1`–`frame_ic_5` files contain complete 3 ps first segments, but their planned 5 ps experiments remain unfinished.

### Li controls completed after the 19 November Release

A read-only server check at **2026-11-22 19:29:27 UTC+08:00** gave the following results:

| Job / case | End time (UTC+08:00) | Technical status | Coverage | S–C mean / minimum / final, Å | Frames with S–C <2 Å |
|---|---|---|---|---:|---:|
| 8572577 / `r1150_original` | 2026-11-20 03:52:44 | COMPLETED, exit 0:0, CHECKS PASS | 501 frames, 0–5000 fs | 3.554832 / 2.718121 / 3.781491 | 0/501 |
| 8572578 / `r2200_original` | 2026-11-20 02:27:27 | COMPLETED, exit 0:0, CHECKS PASS | 501 frames, 0–5000 fs | 3.373292 / 2.518223 / 3.043609 | 0/501 |

Each run has 5000 MD steps, 5001 converged SCFs with no failures, 5001 finite energy rows and a final step-5000 restart. These checks confirm numerical completeness and the distance statistics above. Bond assignment still requires structural review: the current records have `structural_verdict: null` and `human_review: PENDING`.

The verification summary, scheduler receipt and original technical checks are in `evidence/li2025_controls_20261122/`. The 19 November Li Release still contains the older 294-frame and 311-frame prefixes. Keep that Release and its manifest unchanged; publish the complete files as a new snapshot.

## 4. Next work

### WP1 — Release the complete Li controls

No new simulation is needed. Use collection `li2025_hbcontrol_20261118/`, jobs 8572577/8572578 and the hashes in `metadata/HANDOFF_STATUS.json`.

1. Transfer the complete XYZ, energy, stdout, final restart, input, coordinates, run record, `CHECKS.json` and distance-analysis files.
2. Check their hashes before editing any public manifest.
3. Create a new dataset snapshot with 501 frames and steps 0–5000 for both cases.
4. Update `data/DATASETS.json`, `metadata/DATA_FILES.json`, the transfer records and the relevant documentation together. Preserve the old Release and its regression references.
5. Add full-trajectory replay checks and publish a new Release with archive and per-file SHA-256 records.

**Validation:** a clean checkout downloads both complete trajectories, all hashes pass and both replays report 501 frames. The old Release remains available unchanged.

### WP2 — Finish the Li three-arm comparisons

Compare these matched groups:

- `r1150_original`, `r1150_oc_donor`, `r1150_reverse`
- `r2200_original`, `r2200_oc_donor`, `r2200_reverse`

Use full 0–5 ps coverage and the fixed windows 0–0.1, 0.1–1 and 1–5 ps. Include S–C/C–O, carbon geometry, actual water contacts, orientation, signed-GDS depth and the relevant energy records. Save the commands, source hashes and machine-readable outputs; link any S–C-versus-time figure back to those outputs.

The labels describe the initial preparation, not permanent network states. Adjacent frames are correlated, and one trajectory per condition does not establish a general water-orientation effect.

**Deliverables:** complete and fixed-window tables for both groups, with traceable inputs for all six trajectories and a clear separation between observation and interpretation.

### WP3 — Review the structures

Use the working project's three-layer procedure: deterministic geometry screen, model-assisted visual screening and human review. Store the review files and reviewer decision with the trajectory hash. Record any disagreement between the model and the reviewer.

A structural statement is ready for the public report only after this review. Technical PASS and a distance-threshold count are not substitutes for it.

### WP4 — Decide whether to continue Liu

This assessment does not restart sampling.

1. Request the original reactant-complex, input and restart records where possible. Otherwise, document every reconstructed difference, including wall positions and force constants.
2. Inspect component-distance sampling before starting longer metadynamics. A changed CV or wall scheme is a separate experiment with a fresh bias history.
3. If the five frame-start tests are to be completed, validate the step-3000 checkpoints and fix the restart-filename wrapper mismatch before submission. Preserve the 3 ps first segments and identify the short recovery fragments separately.
4. Use independently prepared starts where replication is claimed. A+B and B-only views of the same trajectory count once.

Before any continuation, write down the scientific question, required inputs, coverage target, acceptance criteria and resource cost. A free-energy comparison also needs convergence and method-equivalence checks.

## 5. Reporting and future handoffs

For each work session, record the date, commit, case/job IDs, input and output hashes, and whether the work was local analysis, file transfer or an HPC write. Use a repository issue for data requests and include the collection ID and intended comparison.

Mark trajectories as biased, unbiased or mixed by stage and retain their preparation history. Keep historical outputs, frozen inputs, thresholds and Release assets unchanged. MLP results remain separate from the traditional CP2K evidence. Older Liu orientation and relay outputs with documented binning, sign or periodic-boundary limitations should remain labelled as historical diagnostics.

Update this file or add a dated successor when the status changes. README, progress reports, manifests, citation metadata and Release notes should all describe the same snapshot.

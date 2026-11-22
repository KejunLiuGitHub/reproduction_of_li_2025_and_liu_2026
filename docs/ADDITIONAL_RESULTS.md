# Additional results available on request

This catalogue describes calculations and supporting files not included in full in the public downloads. To request a subset, [open a repository issue](https://github.com/KejunLiuGitHub/reproduction_of_li_2025_and_liu_2026/issues/new) and give the collection ID, case/job or segment, file types, and the comparison you want to make. For a trajectory request, specify the complete run or time interval together with its inputs and preparation history.

Sizes are decimal MB/GB at the stated inventory time. They include large logs and intermediate files where indicated. Collection paths are identifiers within the research archive, not paths expected in a public checkout. Selected analyses and trajectories already available in Releases do not need to be requested again.

## Li 2025

The Li inventory was observed on 19 November 2026 at 21:52:08 UTC+08:00: 2,102 files and 12.679 GB across 16 directories. The public scientific selection contains 583 as-run files plus derived technical traces. The list below counts only files outside that 583-file selection; their recorded sizes precede retrieval of the two active controls.

### LI-01 — New-velocity controls and R-started metaD

- **Collection:** `li2025_branches_20261112/`
- **Not uploaded:** 102 files, 1,234,485,855 bytes (1234.49 MB).

Full logs, wavefunctions and additional restart records for 8569493–8569495: reactant-started metaD, author-product new-velocity control and new-velocity release from the relaxed configuration.

### LI-02 — Duplicate preparation/dry-run package

- **Collection:** `li2025_branches_20261112_dryrun/`
- **Not uploaded:** 14 files, 115,095 bytes (0.12 MB).

Input staging, dry-run validation and duplicate preparation records. They document setup checks but contain no additional production trajectory. Request only for reconstructing the preparation workflow.

### LI-03 — Original bulk/interface runs and early releases

- **Collection:** `li2025_direct/`
- **Not uploaded:** 271 files, 3,371,894,825 bytes (3371.89 MB).

Full CP2K and scheduler logs, wavefunctions, restart backups and supporting files for the first 20 ps bulk/interface chains and their early unbiased releases. The selected release already contains the saved coordinate/energy/CV records needed for the reported comparisons. Request these additional files to inspect SCF histories or reconstruct a particular restart.

### LI-04 — Original interface chain, 20–30 ps

- **Collection:** `li2025_extend_20261111/`
- **Not uploaded:** 35 files, 555,779,719 bytes (555.78 MB).

Additional logs, wavefunctions and restart material for the continuation associated with job 8568711. The selected complete saved trajectory and bias history are already public.

### LI-05 — Original interface chain, 30–40 ps

- **Collection:** `li2025_extend_20261113/`
- **Not uploaded:** 34 files, 555,868,187 bytes (555.87 MB).

Additional logs, wavefunctions and restart material for job 8571154. This interval contains the source of the later 40 ps selected-state release.

### LI-06 — Duplicate preparation/dry-run package

- **Collection:** `li2025_extend_20261113_dryrun/`
- **Not uploaded:** 7 files, 286,911 bytes (0.29 MB).

Input staging, dry-run validation and duplicate preparation records. They document setup checks but contain no additional production trajectory. Request only for reconstructing the preparation workflow.

### LI-07 — Original interface chain, 40–50 ps

- **Collection:** `li2025_extend_20261115/`
- **Not uploaded:** 34 files, 555,872,199 bytes (555.87 MB).

Additional logs, wavefunctions and restart material for job 8572387, including the return to longer C-S distances. The selected release retains the complete saved segment.

### LI-08 — Duplicate preparation/dry-run package

- **Collection:** `li2025_extend_20261115_dryrun/`
- **Not uploaded:** 7 files, 357,669 bytes (0.36 MB).

Input staging, dry-run validation and duplicate preparation records. They document setup checks but contain no additional production trajectory. Request only for reconstructing the preparation workflow.

### LI-09 — Two unrotated original-coordinate controls

- **Collection:** `li2025_hbcontrol_20261118/`
- **Not uploaded:** 54 files, 434,447,868 bytes (434.45 MB).

Additional output/restart material for 8572577/8572578. The `v2026.11.19` public snapshot stops at 2.930/3.100 ps. A later read-only verification found both runs complete at 5 ps, as recorded in [HANDOFF.md](../HANDOFF.md) and `metadata/HANDOFF_STATUS.json`. Until those complete assets appear in a new Release, request this collection with its own retrieval time and coverage and verify the handoff hashes.

### LI-10 — Four local-water orientation tests

- **Collection:** `li2025_hbintervention_20261116/`
- **Not uploaded:** 127 files, 1,295,635,934 bytes (1295.64 MB).

Full output, wavefunctions and segmented restart files for 8572523–8572526. The selected download includes all four complete saved 5 ps trajectories. Request the whole paired group when comparing donor/reverse initial conditions.

### LI-11 — Two selected-state 5 ps releases

- **Collection:** `li2025_hbnetwork_20261116/`
- **Not uploaded:** 66 files, 647,756,258 bytes (647.76 MB).

Full output and restart material for 8572454/8572455. Both starts were already short-C-S configurations drawn from metaD; these records describe subsequent retention, not formation from separated reactants.

### LI-12 — Reactant-started metaD continuation, 10–20 ps

- **Collection:** `li2025_r_extend_20261113/`
- **Not uploaded:** 35 files, 555,818,094 bytes (555.82 MB).

Additional logs and restart material for job 8571156 and its bias history. This is a separate chain from the original 50 ps interface run.

### LI-13 — Duplicate preparation/dry-run package

- **Collection:** `li2025_r_extend_20261113_dryrun/`
- **Not uploaded:** 7 files, 149,664 bytes (0.15 MB).

Input staging, dry-run validation and duplicate preparation records. They document setup checks but contain no additional production trajectory. Request only for reconstructing the preparation workflow.

### LI-14 — Author-snapshot P/R controls

- **Collection:** `li2025_recovery_20261111/`
- **Not uploaded:** 87 files, 761,200,086 bytes (761.20 MB).

Full output and electronic-state/restart files for public-product bulk/interface controls and the 5 ps reactant pre-equilibration (8568510, 8568511 and 8568513). These can resolve detailed numerical and source-preparation questions.

### LI-15 — Restrained preparation and unbiased release

- **Collection:** `li2025_solvation_20261111/`
- **Not uploaded:** 56 files, 490,076,441 bytes (490.08 MB).

Full numerical logs and restart material for job 8568949: 2 ps of S-C restraint followed by 5 ps without that restraint. Request preparation and release together; the preparation history is part of the interpretation.

### LI-16 — Separate static-calculation branch

- **Collection:** `li2025_static/`
- **Not uploaded:** 583 files, 1,851,512,756 bytes (1851.51 MB).

Inputs, optimizations and related output outside the selected dynamics archive. This branch was not used to establish either reproduction result. Request charge/multiplicity, final geometry, termination and frequency records for any proposed stationary-point comparison.

These additional Li files are chiefly full verbose output, wavefunctions and restart backups. They can support SCF inspection or exact restart reconstruction, but are not additional independent simulations. The author workbook and article PDFs should be obtained from their published sources.

## Liu 2026

File sizes below come from a read-only inventory observed on **2026-11-19T22:43:49+08:00**. Its scope was the traditional Liu interface, bulk and proton-scan directories: **5,210 files / 24,984,938,167 bytes (24.985 GB)**. This storage total includes verbose logs, wavefunctions, restart backups, duplicated histories and partial attempts. It is not the size of a recommended scientific download. The file-level inventory excludes W10 and machine-learned-potential collections. A separate top-level W10 storage reading is reported below.

Paths are collection identifiers relative to the research archive; they do not expose a network location or promise that a path exists in the public checkout. An inventory shows a file exists, not that every calculation in the directory reached its intended endpoint. Contents and hashes should be checked when a requested subset is assembled.

### Full coordinate histories of the stopped interface metadynamics chains

The selected release includes the final archived COLVAR/HILLS histories and key inputs. The following coordinate files were left out: **8 files; 240,467,520 bytes (240.47 MB)** in total.

| Collection/file ID | Recorded bytes |
|---|---:|
| `repl_liu2026/iface/hms_iface_resume_seg3/hms_iface-pos-1.xyz` | 31,586,400 |
| `repl_liu2026/iface/hms_iface_resume_seg4/hms_iface-pos-1.xyz` | 31,586,400 |
| `repl_liu2026/iface/hms_iface_resume_seg1/hms_iface-pos-1.xyz` | 31,586,400 |
| `repl_liu2026/iface/hms_iface_resume_seg2/hms_iface-pos-1.xyz` | 29,300,880 |
| `repl_liu2026/iface/metad_seed_step890/hms_iface/hms_iface-pos-1.xyz` | 35,027,520 |
| `repl_liu2026/iface/hms_iface/hms_iface-pos-1.xyz` | 11,658,720 |
| `repl_liu2026/iface/hms_iface_resume_seg5/hms_iface-pos-1.xyz` | 20,158,800 |
| `repl_liu2026/iface/metad_seed_step340/hms_iface/hms_iface-pos-1.xyz` | 49,562,400 |

These cover the original interface chain and its continuation segments, plus the frame-340 and frame-890 branches. The last running jobs were 8554932, 8562431 and 8562419; their recorded times at cancellation were 60.650, 19.290 and 13.610 ps. The eight files are stored segments, not eight independent replicas. Original segments and cumulative CV/bias files can have boundary overlap. A coordinate request should include the matching energy files, segment inputs, source coordinates, HILLS/COLVAR and required checkpoint metadata so that the chronological chain can be reconstructed without inventing missing steps.

These histories may help diagnose solvent-contact loss and which component distances move the five-distance CV. They do not currently support a converged interfacial barrier or a successful reproduction. No new simulation is implied by making the stored files available.

### Source orientation/preparation trajectory

- **Collection:** `repl_liu2026/iface/orient_iface_10ps/`.
- **Inventory:** 58 files; 764,875,434 bytes (764.88 MB) for the entire directory, including logs and restart material.
- **Use:** inspect the source history of the five frame-start tests, its restraints, solvent contacts and choice of 1–5 ps starting configurations.
- **Limit:** the five published B1 trajectories are complete 3 ps segments, but their planned 5 ps tests were not finished. Their starting configurations come from one source history; they are not automatically independent equilibrium samples. The complete source trajectory is not part of the selected release.

### Other unbiased starts and incomplete continuations

- **Collections:** `repl_liu2026/iface/unbiased3_r2/`, `unbiased3_r3/`, `unbiased_sc320b/`, `unbiased_sc320b_r2/`, `unbiased_sc320b_r3/`, and `unbiased_frame_r1/` through `unbiased_frame_r3/` where present.
- **Inventory of those exact directories:** 293 files; 2,130,126,318 bytes (2130.13 MB).
- **Use:** sensitivity to initial RC construction, constrained S–C geometry, selected explicit-solvent frames and velocity initialization.
- **Limit:** many tests cover 3 ps, which is shorter than the later 5 ps observation plan. Their preparation histories differ. They should be compared as named protocols, not pooled into an equilibrium reaction probability.

The small B2 fragments from the five published frame-start cases are also retained outside the selection. Initial jobs 8554943–8554947 ended after a restart-filename assertion even though the B1 calculations and step-3000 checkpoints were present. Recovery jobs 8562388–8562392 were subsequently cancelled; the recorded B2 endpoints of the first three were only steps 3052, 3051 and 3049. Request the individual fragment and its parent checkpoint together if diagnosing continuation behaviour. These fragments are not completed 5 ps results.

### RC construction and interrupted pre-equilibration

- **Collections:** `repl_liu2026/iface/rc_eq*`, `rc_from_eq*`, `rc_sc320*`, and `rc_equil10ps/` where present.
- **Inventory for the listed name families:** 227 files; 911,233,155 bytes (911.23 MB).
- **Contents/use:** restrained solvent preparation, small-cluster optimization/embedding, source-frame searches, input construction records and interrupted pre-equilibration.
- **Limit:** these are preparation steps. A short contact under a restraint or an optimization intermediate is not an unbiased formation event. The 10 ps pre-equilibration attempt 8554980 ended on a wall-clock limit, with energy through step 3225 and a restart at step 3200; the shorter saved checkpoint must be used when discussing recoverability.

### Additional controlled placements, rotations and proton-prepositioning

- **Orientation/placement families:** `repl_liu2026/iface/orientation_placement_*`, `orient_bomd_*`, `manual_*`, `fixed_starts/`.
- **Inventory:** 528 files; 1,615,699,933 bytes (1615.70 MB) across these families; this includes short/partial attempts and source packages.
- **Additional proton-prepositioning families:** `repl_liu2026/iface/proton_preseat_*`; 93 files; 615,905,018 bytes (615.91 MB).
- **Use:** reconstruct which initial orientations, distances, local proton placements and clash checks were tried, and identify why a specific candidate was continued or stopped.
- **Limit:** these are controlled initial-condition searches, not spontaneous-equilibrium ensembles. Some trials are incomplete; names such as “product” or “prebond” describe construction intent and do not establish the chemical identity of a trajectory. Publish or compare the full named batch when drawing a batch-level conclusion.

The four selected D15 A/B arms already provide the complete two-seed retention/control comparison. Their unselected raw segment files, checkpoints and verbose logs can be requested for restart reconstruction. The selected full A+B and phase-B-only trajectories overlap intentionally and must not be counted as independent samples.

### Bulk metadynamics and other bulk comparisons

- **Collections:** `repl_liu2026/bulk/hms_bulk/` and `repl_liu2026/bulk/hmsi_bulk/`.
- **Inventory:** 74 files; 1,669,824,217 bytes (1669.82 MB) for these two directories, including large logs and wavefunctions.
- **Use:** inspect early bulk five-distance sampling and compare the coordinate behaviour with the interface chains.
- **Limit:** these chains were stopped before their planned 100 ps production target. The historical stop record gives approximately 11.98 ps for HMS and 10.67 ps for HMSi, with earlier saved checkpoints. They do not provide the completed bulk free-energy reference needed to validate interfacial acceleration.

Separate bulk guided or unbiased tests are recorded under `repl_liu2026/bulk/`. Their purpose is protocol comparison; an applied force remains a bias even when used in both environments. Request the named protocol and all its replicas rather than treating the whole bulk directory as a single homogeneous data set.

### Proton-coordinate and static geometry scans

- **Collection:** `repl_liu2026/proton_scan/`.
- **Inventory:** 237 files; 23,685,797 bytes (23.69 MB).
- **Use:** examine local proton-coordinate or constrained geometry responses in selected clusters/configurations.
- **Limit:** a constrained static energy scan is not a solvent-averaged free-energy surface. Intermediate optimization energies cannot be substituted for converged stationary-point energies. Normal termination, final geometry, charge/multiplicity and any required frequency checks must be supplied for a particular static comparison.

Other S–C scans and RC geometry calculations appear under `repl_liu2026/iface/sc_scan*`. They are separate diagnostics and do not establish the reported dynamical barrier.

### W10 cluster search and higher-level quantum chemistry

- **Collections:** `repl_liu2026/w10/` and the local audit collection `results/repl_liu2026/w10_audit/`.
- **Remote data size:** 31,461,567,471 bytes (31.462 GB) in a separate top-level storage observation on 19 November; file-level counts and contents were not inventoried. The local audit directory contains 1376 files totaling 4,628,847 bytes; this smaller number describes audit records only.
- **Contents/use:** cluster candidate generation, optimization funnels, recovery records and proposed higher-level energy comparisons for HMS/HMSi.
- **Limit:** the project used a candidate-search funnel with GFN2-xTB screening and SMD-containing optimization stages; this differs from the paper's static W10 pathway and parts of its search procedure. It cannot reproduce the paper's 4.5 kcal/mol product-energy difference by matching an unrelated funnel ranking. No complete W10/DLPNO reproduction result is claimed in the public progress report.

### Machine-learned-potential exploration

- **Known collection IDs:** `repl_liu2026/mlp_bench/` and `repl_liu2026/mlp_data/`.
- **Size and complete output inventory:** **not inventoried** here.
- **Use:** independent model benchmarking, labelled-configuration generation and exploratory sampling.
- **Limit:** these are outside the traditional CP2K/PBE-D3 reproduction evidence. MLP trajectories cannot substitute for the omitted CP2K production, and their energies or free energies must not be merged with the CP2K or static W10 results. A request should specify model/checkpoint, training-data provenance and reference electronic-structure labels; novel reactive configurations need their own reference calculations.

### Historical diagnostics and large technical files

Historical `interface_orientation` and `proton_relay` outputs are retained in `results/repl_liu2026/`. They were not selected as new quantitative mechanism evidence: the older orientation script uses an unsmoothed 1 Å binning convention and has a lower-interface sign inconsistency, and the older proton-event script uses distances without periodic boundaries. The public reanalysis instead uses the retained core periodic-distance and interface routines. Old files should remain identifiable as outputs of their original algorithms.

Within the inventoried traditional scope, filenames ending in `.out` or `.log` account for **367 files; 18,982,902,429 bytes (18982.90 MB)**; filenames containing `.wfn` account for **410 files; 4,665,120,204 bytes (4665.12 MB)**. These categories overlap with the collections above and must not be added to their totals. Full stdout can be useful for detailed SCF diagnosis, and wavefunctions/checkpoints for restarting a specific calculation. They are omitted from the selected release because they are large and are not required for basic trajectory reanalysis. Transport warnings, credentials and unrelated personal material are not part of a scientific request package.

When requesting additional files, specify the collection ID, desired time range or whole run, the observable to be checked, and whether restart capability is needed. The response should include the exact file list, sizes, hashes, preparation/bias history and any unfinished segments. Availability of a stored file does not by itself resolve the remaining scientific comparison.

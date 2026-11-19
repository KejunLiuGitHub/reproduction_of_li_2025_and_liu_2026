# Li 2025 data layout and selection

The data release is a selected research snapshot from 19 November 2026. It is organized by computational role, with complete saved trajectories for the included runs. Selection was not based on which frames had short C–S distances.

The archive contains 583 as-run files and 74 technical-trace/index files, 657 files in total. After extraction:

```text
data/server/li2025_direct/                      original bulk/interface and early releases
data/server/li2025_extend_20261111/             interface 20–30 ps
data/server/li2025_extend_20261113/             interface 30–40 ps
data/server/li2025_extend_20261115/             interface 40–50 ps
data/server/li2025_recovery_20261111/           public P/R controls
data/server/li2025_solvation_20261111/          restrained preparation and release
data/server/li2025_branches_20261112/           new-velocity controls and R-started metaD
data/server/li2025_r_extend_20261113/           R-started metaD 10–20 ps
data/server/li2025_hbnetwork_20261116/           two selected-state 5 ps releases
data/server/li2025_hbintervention_20261116/      four water-orientation experiments
data/server/li2025_hbcontrol_20261118/           two unrotated controls
data/technical_traces/                         compact records extracted from CP2K stdout
```

`data/DATASETS.json` lists 46 XYZ trajectories, including short smoke runs. Each entry gives the matching energy and input files, box, bias flag and saved-frame coverage. Continuation segments retain their original step numbers; the original and R-started metadynamics chains have separate bias histories.

The original controls are live snapshots captured on 19 November 2026, with retrieval recorded at 21:58:59 UTC+08:00: R1150 contains 294 complete frames through 2930 fs; R2200 contains 311 through 3100 fs. Their energy, restart and stdout records may end at different steps. They are not 5 ps results. Completed-file hashes were compared before and after retrieval; active files are identified separately in `metadata/TRANSFER_MANIFEST.json`.

## Included

- Scientific inputs and source coordinates, seeds, scripts and parameter manifests.
- Complete saved XYZ and energy sequences for the selected runs.
- HILLS/COLVAR and final restart files needed to inspect bias history and continuation state.
- Existing numerical summaries and run records, including unsuccessful or short runs.
- Compact CP2K traces with header/tail information and SCF/MD markers.

## Omitted

Wavefunctions, redundant restart backups, duplicate dry-runs, operating-system caches, the separate unregistered static branch and full verbose stdout are omitted. The traces are explicitly labelled **EXCERPT_NOT_FULL_LOG** and carry original-log or captured-prefix hashes. SSH connection warnings have been removed from these derived traces; the original-log and captured-prefix hashes are preserved. They support inspection of the recorded run but do not replace full stdout for every historical diagnostic. Full verbose logs remain in the working archive.

The author workbook and paper PDFs are obtained from their original sources. No email, credentials, unrelated project data or private collaboration notes are included.

## Integrity and size

`metadata/DATA_FILES.json` lists each released file's size and SHA-256. `metadata/RELEASE_ASSETS.json` lists the compressed archive checksum. Use `python tools/fetch_data.py --verify-only` to check a local copy. `tools/fetch_data.py --archive /path/to/archive.tar.gz` also supports an offline copy.

Data are distributed as a GitHub Release attachment to keep the Git repository small. GitHub's [release documentation](https://docs.github.com/en/repositories/releasing-projects-on-github/about-releases) describes the attachment mechanism.

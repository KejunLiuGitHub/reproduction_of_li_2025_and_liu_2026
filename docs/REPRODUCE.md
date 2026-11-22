# Reproduce the analyses and prepare a fresh run

Use Python 3.12. The replay examples were tested with Python 3.12.10 and NumPy 2.2.6.

```bash
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
```

Download and verify the selected data with `tools/fetch_data.py`. Cases and their paths are
listed in `data/DATASETS.json`. A case's `status` describes the captured data. In Release
`v2026.11.19`, the original-coordinate controls are live prefixes at the stated time. Their later
completion is recorded in [HANDOFF.md](../HANDOFF.md), but the complete assets require the new
Release described there as WP1.

## Replay a case

```bash
python tools/reproduce_case.py --case natural_original40 --out analysis/natural_original40
python tools/reproduce_case.py --case r1150_original --out analysis/r1150_original
python tools/reproduce_case.py --case p_bulk_smoke --out analysis/p_bulk_smoke
```

The command runs the original distance, solvation and network programs. It reads the cell and
bias history from the case manifest, selects the bulk or interface solvation mode, and retains
every complete XYZ frame. Network analysis is enabled by default; use `--no-network` to omit it.
The output directory must be new. Input data are not modified.

For a fixed common comparison window, add `--through-step 2570` (or the stated window endpoint).
`--through-step 100` retains 11 frames at 0,10,…,100 fs. The default uses the full available
trajectory. Selection writes separate copies and records both original and selected hashes,
requested cutoff and actual frame/energy coverage.

`SUMMARY.json` reports SC/CO distances, temperature, water coordination, hydrogen-bond counts
at three thresholds, and interface signed-GDS/normal orientation when applicable. Bulk summaries
have no GDS surface. `COMMANDS.json`, `SOURCE_SHA256.json` and `OUTPUT_SHA256.json` record the
commands and hashes. Original analyzer outputs, complete-frame copies and any unfinished tail
are also retained. A missing historical job number stays missing; its network label is `archive:<case>`.
`DISTANCE_ANALYSIS.json` uses the case's actual bias history and neutral distance-span labels;
`DISTANCE_ANALYSIS_RAW.json` preserves the unchanged legacy output.

XYZ and energy endpoints are reported separately. Replaying a live prefix does not complete the
underlying simulation. Compare multiple arms only over common coverage; the raw distance recorder's
legacy “dissociation” field is an SC-distance span, and the HB spans are sampled observations.
Preparation and bias history are included in the public summary.

## Prepare the two original-coordinate controls

Install CP2K 9.1 and its data files separately. Point to a directory containing `BASIS_MOLOPT`,
`GTH_POTENTIALS` and `dftd3.dat`. Use a path without whitespace.

```bash
python tools/prepare_unbiased_run.py --case r1150_original \
  --cp2k-data-dir /path/to/cp2k-9.1/data --out new_runs/r1150_original
python tools/prepare_unbiased_run.py --case r2200_original \
  --cp2k-data-dir /path/to/cp2k-9.1/data --out new_runs/r2200_original
```

The tool checks the selected input and coordinates against the frozen package. It copies the
coordinates exactly and changes only the three CP2K data-file paths and the output `PROJECT`
name. Seeds, 5000 steps, cell, electronic method, thermostat and absence of bias remain unchanged.
`PREPARATION.json` records old/new hashes and the preparation history; `INPUT_PATH_CHANGES.diff`
shows every edited line. No calculation starts during preparation, and no old run certificate is copied.
The provided CP2K data files are hashed; use the recorded version when comparing numerical results.

Configure your executable and MPI launcher, then use the shell example from the repository root:

```bash
export CP2K_BINARY=/path/to/cp2k.psmp
export MPI_LAUNCHER=mpirun MPI_TASKS=64
bash examples/run_prepared.sh new_runs/r1150_original
```

The example uses 64 MPI tasks and one OpenMP thread, records engine version and return code,
and refuses to launch the same prepared directory twice. `examples/run_prepared.slurm` supplies
a scheduler template; set your own account, partition and environment. Neither example connects
to a particular cluster. New calculations produce their own logs and numerical records.

## Numerical checks for a new calculation

For a complete 5 ps run, inspect 5000 consecutive MD steps, 5001 converged SCFs with zero failed
SCFs, 501 complete 296-atom XYZ frames, 5001 finite energy rows, normal program termination,
and a final restart containing 296 positions and velocities at step 5000. Check first-frame
coordinates against the prepared source, energy comments against the energy table, and preserve
all new input/output hashes. These checks concern the run's numerical completeness; structural
interpretation uses the separately documented geometry analysis.

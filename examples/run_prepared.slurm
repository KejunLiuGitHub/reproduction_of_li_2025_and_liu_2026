#!/usr/bin/env bash
#SBATCH --job-name=li296-original
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=48:00:00
# Add your site's account/partition settings when submitting this example.
set -euo pipefail
: "${PUBLIC_REPO:?Set PUBLIC_REPO to this repository's absolute path}"
: "${PREPARED_DIR:?Set PREPARED_DIR to the new prepared directory}"
: "${CP2K_BINARY:?Set CP2K_BINARY to the CP2K 9.1 executable}"
export MPI_TASKS=64 OMP_NUM_THREADS=1
# Use your site's supported MPI launcher; e.g. MPI_LAUNCHER=srun MPI_NP_FLAG=-n.
bash "$PUBLIC_REPO/examples/run_prepared.sh" "$PREPARED_DIR"

#!/bin/bash
#SBATCH --job-name=li25-recovery
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=32G
#SBATCH --time=24:00:00
set -euo pipefail
BASE=/public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/li2025_recovery_20261111
CASE=${1:?case required}
export OMP_NUM_THREADS=1
export I_MPI_ROOT=/public/soft/intel/oneapi/mpi/latest
export LD_LIBRARY_PATH=/public/home/temp/apps/plumed-2.7.3-shared/lib:/public/soft/intel/oneapi/mkl/latest/lib/intel64:/public/soft/intel/oneapi/compiler/latest/lib:/public/soft/intel/oneapi/tbb/latest/lib:${LD_LIBRARY_PATH:-}
/public/home/temp/miniconda3/bin/python3 "$BASE/run_pilot.py" --manifest pilot_manifest.json --case "$CASE" --ntasks "$SLURM_NTASKS"
/public/home/temp/miniconda3/bin/python3 "$BASE/analyze_pilot.py" "$CASE"

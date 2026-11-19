#!/bin/bash
#SBATCH --job-name=li25-branches
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=32G
#SBATCH --time=48:00:00
set -euo pipefail
BASE=/public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/li2025_branches_20261112
CASE=${1:?branch case required}
export OMP_NUM_THREADS=1
export I_MPI_ROOT=/public/soft/intel/oneapi/mpi/latest
export LD_LIBRARY_PATH=/public/home/temp/apps/plumed-2.7.3-shared/lib:/public/soft/intel/oneapi/mkl/latest/lib/intel64:/public/soft/intel/oneapi/compiler/latest/lib:/public/soft/intel/oneapi/tbb/latest/lib:${LD_LIBRARY_PATH:-}
/public/home/temp/miniconda3/bin/python3 "$BASE/run_branch.py" --case "$CASE" --ntasks "$SLURM_NTASKS"

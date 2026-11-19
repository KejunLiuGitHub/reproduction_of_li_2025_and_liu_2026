#!/bin/bash
#SBATCH --job-name=li25-direct-bulk-r1fix1
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/logs/li25-direct-bulk-r1_%j.out
#SBATCH --error=/public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/logs/li25-direct-bulk-r1_%j.err

# Li 2025 Figure 6b R_bulk -> pathway-3 direct C-S addition, first production replica.
# 64 MPI ranks x 1 thread, no GPUs. Environment is baked in here because hpc sbatch
# does not forward the submitting shell's environment.

set -euo pipefail
export LD_LIBRARY_PATH=/public/home/temp/apps/plumed-2.7.3-shared/lib:/public/soft/intel/oneapi/mkl/latest/lib/intel64:/public/soft/intel/oneapi/compiler/latest/lib:/public/soft/intel/oneapi/tbb/latest/lib:${LD_LIBRARY_PATH:-}
export I_MPI_ROOT=/public/soft/intel/oneapi/mpi/latest
export OMP_NUM_THREADS=1

BASE=/public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/li2025_direct
WORK=$BASE/bulk_r1_fix1
mkdir -p "$WORK" /public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/logs
cd "$WORK"

if [ -s "$WORK/run_record.json" ] && grep -q '"status": "COMPLETED"' "$WORK/run_record.json"; then
    echo "SKIP: $WORK already COMPLETED (remove run_record.json to force a rerun)"
    exit 0
fi

python3 "$BASE/run_cp2k.py" --case bulk_r1 --workdir "$WORK" --ntasks "$SLURM_NTASKS"

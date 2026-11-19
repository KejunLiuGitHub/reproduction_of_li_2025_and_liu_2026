#!/bin/bash
#SBATCH --job-name=li25-direct-bulk-seg2
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/logs/li25-direct-bulk-seg2_%j.out
#SBATCH --error=/public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/logs/li25-direct-bulk-seg2_%j.err

# second 10 ps segment of bulk_r1 (parent job 8562691, workdir bulk_r1_fix1).
# The parent checkpoint, XYZ, inputs, HILLS and COLVAR are copied by the runner
# into bulk_r1_seg2 under non-colliding names and SHA-checked before launch.
# Environment is baked in here because hpc sbatch does not forward the shell env.

set -euo pipefail
export LD_LIBRARY_PATH=/public/home/temp/apps/plumed-2.7.3-shared/lib:/public/soft/intel/oneapi/mkl/latest/lib/intel64:/public/soft/intel/oneapi/compiler/latest/lib:/public/soft/intel/oneapi/tbb/latest/lib:${LD_LIBRARY_PATH:-}
export I_MPI_ROOT=/public/soft/intel/oneapi/mpi/latest
export OMP_NUM_THREADS=1

BASE=/public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/li2025_direct
WORK=$BASE/bulk_r1_seg2
PARENT=$BASE/bulk_r1_fix1
mkdir -p "$WORK" /public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/logs
cd "$WORK"

[ -s "$BASE/seg2_manifest.json" ] || { echo "missing $BASE/seg2_manifest.json" >&2; exit 1; }
[ -s "$PARENT/bulk_r1-1.restart" ] || { echo "missing parent restart in $PARENT" >&2; exit 1; }
[ -s "$PARENT/run_record.json" ] || { echo "missing parent run_record.json" >&2; exit 1; }
grep -q '"status": "COMPLETED"' "$PARENT/run_record.json" \
    || { echo "parent is not COMPLETED: $PARENT/run_record.json" >&2; exit 1; }

if [ -s "$WORK/run_record.json" ] && grep -q '"status": "COMPLETED"' "$WORK/run_record.json"; then
    echo "SKIP: $WORK already COMPLETED (remove run_record.json to force a rerun)"
    exit 0
fi

python3 "$BASE/run_cp2k.py" --manifest seg2_manifest.json --case bulk_r1_seg2 \
    --workdir "$WORK" --parent-dir "$PARENT" --ntasks "$SLURM_NTASKS"

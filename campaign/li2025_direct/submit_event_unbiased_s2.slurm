#!/bin/bash
#SBATCH --job-name=li25-uev-s2
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/logs/li25-uev-s2_%j.out
#SBATCH --error=/public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/logs/li25-uev-s2_%j.err

# Fully unbiased 3 ps replica 2 from the Li 2025 bulk event frame (S-C = 1.778 A).
# No PLUMED, no EXT_RESTART, no inherited velocities: CP2K draws fresh
# Maxwell-Boltzmann velocities at 298.15 K from &GLOBAL SEED 20261110.
# The three replicas are independent and may be submitted as slots become available.

set -euo pipefail
export LD_LIBRARY_PATH=/public/home/temp/apps/plumed-2.7.3-shared/lib:/public/soft/intel/oneapi/mkl/latest/lib/intel64:/public/soft/intel/oneapi/compiler/latest/lib:/public/soft/intel/oneapi/tbb/latest/lib:${LD_LIBRARY_PATH:-}
export I_MPI_ROOT=/public/soft/intel/oneapi/mpi/latest
export OMP_NUM_THREADS=1

BASE=/public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/li2025_direct
WORK=$BASE/event_unbiased_s2
mkdir -p "$WORK" /public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/logs
cd "$WORK"

[ -s "$BASE/event_unbiased_manifest.json" ] || { echo "missing event_unbiased_manifest.json" >&2; exit 1; }

if [ -s "$WORK/run_record.json" ] && grep -q '"status": "COMPLETED"' "$WORK/run_record.json"; then
    echo "SKIP: $WORK already COMPLETED (remove run_record.json to force a rerun)"
    exit 0
fi

python3 "$BASE/run_cp2k.py" --manifest event_unbiased_manifest.json --case event_unbiased_s2 \
    --workdir "$WORK" --ntasks "$SLURM_NTASKS"

if [ -s "$WORK/bulk_r1_unbiased-pos-1.xyz" ]; then
    python3 "$BASE/analyze_event_unbiased.py" \
        --trajectory "$WORK/bulk_r1_unbiased-pos-1.xyz" \
        --energy "$WORK/bulk_r1_unbiased-1.ener" \
        --label event_unbiased_s2 --box 14.2296,14.2296,14.2296 \
        --out "$WORK/EVENT_UNBIASED_2.json" \
        || echo "WARNING: analyzer failed for event_unbiased_s2" >&2
fi

#!/bin/bash
#SBATCH --job-name=li25-iev-s1
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/logs/li25-iev-s1_%j.out
#SBATCH --error=/public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/logs/li25-iev-s1_%j.err

# Fully unbiased 3 ps replica 1 from the Li 2025 INTERFACE short-contact frame
# (step 10410, S-C = 1.808 A, 15 x 15 x 30 A cell), harvested from the biased
# interface_r1 metaD run. No PLUMED, no EXT_RESTART, no inherited velocities:
# CP2K draws fresh Maxwell-Boltzmann velocities at 298.15 K from &GLOBAL SEED 20261112.
# The analyzer must be given the interface cell 15,15,30 - never the bulk 14.2296 cube.

set -euo pipefail
export LD_LIBRARY_PATH=/public/home/temp/apps/plumed-2.7.3-shared/lib:/public/soft/intel/oneapi/mkl/latest/lib/intel64:/public/soft/intel/oneapi/compiler/latest/lib:/public/soft/intel/oneapi/tbb/latest/lib:${LD_LIBRARY_PATH:-}
export I_MPI_ROOT=/public/soft/intel/oneapi/mpi/latest
export OMP_NUM_THREADS=1

BASE=/public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/li2025_direct
WORK=$BASE/interface_event_unbiased_s1
mkdir -p "$WORK" /public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/logs
cd "$WORK"

[ -s "$BASE/interface_event_unbiased_manifest.json" ] || { echo "missing interface_event_unbiased_manifest.json" >&2; exit 1; }

if [ -s "$WORK/run_record.json" ] && grep -q '"status": "COMPLETED"' "$WORK/run_record.json"; then
    echo "SKIP: $WORK already COMPLETED (remove run_record.json to force a rerun)"
    exit 0
fi

python3 "$BASE/run_cp2k.py" --manifest interface_event_unbiased_manifest.json \
    --case interface_event_unbiased_s1 --workdir "$WORK" --ntasks "$SLURM_NTASKS"

if [ -s "$WORK/interface_r1_unbiased-pos-1.xyz" ]; then
    python3 "$BASE/analyze_event_unbiased.py" \
        --trajectory "$WORK/interface_r1_unbiased-pos-1.xyz" \
        --energy "$WORK/interface_r1_unbiased-1.ener" \
        --label interface_event_unbiased_s1 --box 15.0,15.0,30.0 \
        --out "$WORK/INTERFACE_EVENT_UNBIASED_1.json" \
        || echo "WARNING: analyzer failed for interface_event_unbiased_s1" >&2
fi

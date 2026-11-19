#!/bin/bash
#SBATCH --job-name=li25-r-10to20
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=32G
#SBATCH --time=48:00:00
set -euo pipefail
BASE=/public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/li2025_r_extend_20261113
PARENT=/public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/li2025_branches_20261112/r_eq_metad
/public/home/temp/miniconda3/bin/python3 "$BASE/run_continuation.py" --manifest continuation_manifest.json --case r_eq_metad_seg2 --parent-dir "$PARENT" --ntasks "$SLURM_NTASKS"

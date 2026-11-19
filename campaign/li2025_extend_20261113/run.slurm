#!/bin/bash
#SBATCH --job-name=li25-iface-30to40
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=32G
#SBATCH --time=48:00:00
set -euo pipefail
BASE=/public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/li2025_extend_20261113
PARENT=/public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/li2025_extend_20261111/interface_r1_seg3
/public/home/temp/miniconda3/bin/python3 "$BASE/run_continuation.py" --manifest continuation_manifest.json --case interface_r1_seg4 --parent-dir "$PARENT" --ntasks "$SLURM_NTASKS"

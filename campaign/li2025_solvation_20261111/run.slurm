#!/bin/bash
#SBATCH --job-name=li25-solv-relax
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=32G
#SBATCH --time=48:00:00
set -euo pipefail
BASE=/public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/li2025_solvation_20261111
/public/home/temp/miniconda3/bin/python3 "$BASE/run.py" --ntasks "$SLURM_NTASKS"

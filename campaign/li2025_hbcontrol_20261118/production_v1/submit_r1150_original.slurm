#!/bin/bash
#SBATCH --job-name=li25-r1150_original
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=32G
#SBATCH --time=48:00:00
set -euo pipefail
BASE="${LI_WORKFLOW_BASE:?Use the frozen launch wrapper with an absolute package directory}"
/public/home/temp/miniconda3/bin/python3 "$BASE/run.py" --base "$BASE" --case r1150_original --ntasks "$SLURM_NTASKS"

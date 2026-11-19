#!/bin/bash
#SBATCH --job-name=li25-r2200_oc_donor
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=32G
#SBATCH --time=48:00:00
set -euo pipefail
BASE="${LI_WORKFLOW_BASE:?Set absolute uploaded package directory in registered wrapper}"
/public/home/temp/miniconda3/bin/python3 "$BASE/run.py" --base "$BASE" --case r2200_oc_donor --ntasks "$SLURM_NTASKS"

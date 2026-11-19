#!/bin/bash
#SBATCH --job-name=li25-natural_original40
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=32G
#SBATCH --time=48:00:00
set -euo pipefail
BASE="${LI_WORKFLOW_BASE:?Set the absolute uploaded package directory in a registered submission wrapper}"
/public/home/temp/miniconda3/bin/python3 "$BASE/run.py" --base "$BASE" --case natural_original40 --ntasks "$SLURM_NTASKS"

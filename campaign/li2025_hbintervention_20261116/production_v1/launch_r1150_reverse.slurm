#!/bin/bash
#SBATCH --job-name=li25-r1150_reverse
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=32G
#SBATCH --time=48:00:00
set -euo pipefail
export LI_WORKFLOW_BASE="/public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/li2025_hbintervention_20261116/production_v1"
source "$LI_WORKFLOW_BASE/submit_r1150_reverse.slurm"

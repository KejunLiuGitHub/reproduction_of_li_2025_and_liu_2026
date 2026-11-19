#!/bin/bash
#SBATCH --job-name=li25-natural_r19830
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --chdir=/public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/li2025_hbnetwork_20261116/production_v1
#SBATCH --output=/public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/li2025_hbnetwork_20261116/production_v1/launch-%j.out
set -euo pipefail
export LI_WORKFLOW_BASE='/public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/li2025_hbnetwork_20261116/production_v1'
export PYTHONDONTWRITEBYTECODE=1
source "$LI_WORKFLOW_BASE/submit_natural_r19830.slurm"

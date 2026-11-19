#!/bin/bash
#SBATCH --job-name=li25-r_eq_metad
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/li2025_branches_20261112/r_eq_metad_%j.out
#SBATCH --error=/public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/li2025_branches_20261112/r_eq_metad_%j.err
# session: li2025-branches-20261112
source /public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/li2025_branches_20261112/run.slurm r_eq_metad

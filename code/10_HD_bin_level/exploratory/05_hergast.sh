#!/bin/bash
#SBATCH -p caracol
#SBATCH --gpus=1
#SBATCH --mem=100G
#SBATCH --job-name=05_hergast
#SBATCH -c 1
#SBATCH -t 2-00:00:00
#SBATCH -o ../../processed-data/10_HD_bin_level/logs/05_hergast.txt
#SBATCH -e ../../processed-data/10_HD_bin_level/logs/05_hergast.txt

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load hergast/0.0.1

## List current modules for reproducibility
module list

python 05_hergast.py

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/

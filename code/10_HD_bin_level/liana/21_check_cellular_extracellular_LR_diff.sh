#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=50G
#SBATCH --job-name=21_check_cellular_extracellular_LR_diff
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../../processed-data/10_HD_bin_level/no_secondary/liana/logs/21_check_cellular_extracellular_LR_diff.txt
#SBATCH -e ../../../processed-data/10_HD_bin_level/no_secondary/liana/logs/21_check_cellular_extracellular_LR_diff.txt

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load liana_plus

## List current modules for reproducibility
module list

python3 21_5_check_cellular_extracellular_LR_diff_cell_type2.py
python3 21_check_cellular_extracellular_LR_diff.py

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/

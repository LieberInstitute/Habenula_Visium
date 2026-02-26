#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=50G
#SBATCH --job-name=22_heatmap
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../../processed-data/10_HD_bin_level/no_secondary/liana/logs/22_heatmap.txt
#SBATCH -e ../../../processed-data/10_HD_bin_level/no_secondary/liana/logs/22_heatmap.txt

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

python3 22_heatmap_LR_cell_type.py
python3 22_heatmap_LR_cell_type_specific.py

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/

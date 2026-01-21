#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=32G
#SBATCH --job-name=05_secondary_prob
#SBATCH -c 1
#SBATCH -t 8:00:00
#SBATCH -o ../../../processed-data/09_HD_cell_level/new_samples2/refine_segmentations/logs/05_secondary_prob_%a.txt
#SBATCH -e ../../../processed-data/09_HD_cell_level/new_samples2/refine_segmentations/logs/05_secondary_prob_%a.txt
#SBATCH --array=30

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load visium_hd/1.0

## List current modules for reproducibility
module list

python 05_secondary_prob.py

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/

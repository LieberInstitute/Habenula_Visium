#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=64G
#SBATCH --job-name=18_secondary_segmentations
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o ../../processed-data/09_HD_cell_level/new_samples2/logs/18_secondary_segmentations.txt
#SBATCH -e ../../processed-data/09_HD_cell_level/new_samples2/logs/18_secondary_segmentations.txt

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

python 18_secondary_segmentations.py

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/

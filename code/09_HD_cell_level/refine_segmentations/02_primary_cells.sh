#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=32G
#SBATCH --job-name=02_primary_cells
#SBATCH -c 1
#SBATCH -t 8:00:00
#SBATCH -o ../../../processed-data/09_HD_cell_level/probe_fix/refine_segmentations/logs/02_primary_cells_%a.txt
#SBATCH -e ../../../processed-data/09_HD_cell_level/probe_fix/refine_segmentations/logs/02_primary_cells_%a.txt
#SBATCH --array=1-5%5

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

python 02_primary_cells.py

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/

#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=20G
#SBATCH --job-name=13_singleR
#SBATCH -c 1
#SBATCH -t 2-0:00:00
#SBATCH -o ../../processed-data/09_HD_cell_level/logs/13_singleR_%a.txt
#SBATCH -e ../../processed-data/09_HD_cell_level/logs/13_singleR_%a.txt
#SBATCH --array=1-2%2

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load the R module
module load conda_R/4.4.x

## List current modules for reproducibility
module list

Rscript 13_singleR.R

echo "**** Job ends ****"
date

#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=8G
#SBATCH --job-name=15_crawdad_run
#SBATCH -c 2
#SBATCH -t 1:00:00
#SBATCH -o ../../processed-data/09_HD_cell_level/no_secondary/logs/15_crawdad_run_%a.txt
#SBATCH -e ../../processed-data/09_HD_cell_level/no_secondary/logs/15_crawdad_run_%a.txt
#SBATCH --array=9
#SBATCH --exclude=compute-169

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
module load conda_R/4.5

## List current modules for reproducibility
module list

Rscript 15_crawdad_run.R

echo "**** Job ends ****"
date

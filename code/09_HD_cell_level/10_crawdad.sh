#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=10_crawdad
#SBATCH -c 4
#SBATCH -t 2-0:00:00
#SBATCH -o ../../processed-data/09_HD_cell_level/new_samples/logs/10_crawdad_%a.txt
#SBATCH -e ../../processed-data/09_HD_cell_level/new_samples/logs/10_crawdad_%a.txt
#SBATCH --array=1-3%3

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

Rscript 10_crawdad.R

echo "**** Job ends ****"
date

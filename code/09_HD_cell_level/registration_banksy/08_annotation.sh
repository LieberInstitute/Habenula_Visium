#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=16G
#SBATCH --job-name=08_annotation
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../processed-data/09_HD_cell_level/new_samples2/registration_banksy/logs/08_annotation.txt
#SBATCH -e ../../processed-data/09_HD_cell_level/new_samples2/registration_banksy/logs/08_annotation.txt

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

Rscript 08_annotation.R

echo "**** Job ends ****"
date

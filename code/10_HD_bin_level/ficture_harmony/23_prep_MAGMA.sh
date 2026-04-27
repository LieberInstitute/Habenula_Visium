#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=40G
#SBATCH --job-name=23_prep_MAGMA
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../../processed-data/10_HD_bin_level/new_samples2/ficture_harmony/logs/23_prep_MAGMA.txt
#SBATCH -e ../../../processed-data/10_HD_bin_level/new_samples2/ficture_harmony/logs/23_prep_MAGMA.txt

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

Rscript 23_prep_MAGMA.R

echo "**** Job ends ****"
date

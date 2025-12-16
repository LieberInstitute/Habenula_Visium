#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=15G
#SBATCH --job-name=18_LIANA+_NMF
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../../processed-data/10_HD_bin_level/new_samples2/liana/logs/18_LIANA+_NMF_%a.txt
#SBATCH -e ../../../processed-data/10_HD_bin_level/new_samples2/liana/logs/18_LIANA+_NMF_%a.txt
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
module load liana_plus

## List current modules for reproducibility
module list

python3 18_LIANA+_NMF.py

echo "**** Job ends ****"
date
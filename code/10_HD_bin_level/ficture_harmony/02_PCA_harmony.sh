#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=100G
#SBATCH --job-name=02_PCA_harmony
#SBATCH -c 10
#SBATCH -t 2-0:00:00
#SBATCH -o ../../../processed-data/10_HD_bin_level/ficture_harmony/logs/02_PCA_harmony.txt
#SBATCH -e ../../../processed-data/10_HD_bin_level/ficture_harmony/logs/02_PCA_harmony.txt

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

Rscript 02_PCA_harmony.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/

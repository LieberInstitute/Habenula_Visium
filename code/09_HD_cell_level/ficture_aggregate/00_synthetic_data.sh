#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=40G
#SBATCH --job-name=00_synthetic_data
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../../processed-data/09_HD_cell_level/ficture_aggregate/logs/00_synthetic_data.txt
#SBATCH -e ../../../processed-data/09_HD_cell_level/ficture_aggregate/logs/00_synthetic_data.txt

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

Rscript 00_synthetic_data.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/

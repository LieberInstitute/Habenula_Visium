#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=300G
#SBATCH --job-name=02_build_spe_raw
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../processed-data/09_HD_cell_level/new_samples/logs/02_build_spe_raw.txt
#SBATCH -e ../../processed-data/09_HD_cell_level/new_samples/logs/02_build_spe_raw.txt

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

Rscript 02_build_spe_raw.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/

#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=50G
#SBATCH --job-name=06_add_to_spe
#SBATCH -c 1
#SBATCH -t 2:00:00
#SBATCH -o ../../../processed-data/10_HD_bin_level/logs/all_sample_002um/06_add_to_spe.txt
#SBATCH -e ../../../processed-data/10_HD_bin_level/logs/all_sample_002um/06_add_to_spe.txt

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
module load conda_R/4.4

## List current modules for reproducibility
module list

Rscript 06_add_to_spe.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/

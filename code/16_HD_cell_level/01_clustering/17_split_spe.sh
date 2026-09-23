#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=32G
#SBATCH --job-name=17_split_spe
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH -o ../../processed-data/09_HD_cell_level/no_secondary/logs/17_split_spe.txt
#SBATCH -e ../../processed-data/09_HD_cell_level/no_secondary/logs/17_split_spe.txt

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

Rscript 17_split_spe.R

echo "**** Job ends ****"
date

#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=22_top_LR_across_cells
#SBATCH -t 1-00:00:00
#SBATCH -o ../../../processed-data/10_HD_bin_level/no_secondary/liana2/logs/22_top_LR_across_cells.txt
#SBATCH -e ../../../processed-data/10_HD_bin_level/no_secondary/liana2/logs/22_top_LR_across_cells.txt

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

## Edit with your job command
Rscript 22_top_LR_across_cells.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.4.0
## available from http://research.libd.org/slurmjobs/

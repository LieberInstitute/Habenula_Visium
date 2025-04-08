#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=240G
#SBATCH --job-name=05_banksy_embedding
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../processed-data/09_HD_cell_level/probe_fix/logs/05_banksy_embedding_%a.txt
#SBATCH -e ../../processed-data/09_HD_cell_level/probe_fix/logs/05_banksy_embedding_%a.txt
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
module load conda_R/4.4.x

## List current modules for reproducibility
module list

Rscript 05_banksy_embedding.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/

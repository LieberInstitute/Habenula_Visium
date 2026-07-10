#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=20G
#SBATCH --job-name=05_prep_visium_shiny
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../processed-data/12_apps_and_sharing/logs/05_prep_visium_shiny.txt
#SBATCH -e ../../processed-data/12_apps_and_sharing/logs/05_prep_visium_shiny.txt

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

Rscript 05_prep_visium_shiny.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/

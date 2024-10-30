#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=08_plot_SVGs
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH -o ../../processed-data/10_HD_bin_level/logs/08_plot_SVGs.txt
#SBATCH -e ../../processed-data/10_HD_bin_level/logs/08_plot_SVGs.txt

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

Rscript 08_plot_SVGs.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/

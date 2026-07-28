#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=02_cor_heatmap
#SBATCH -c 1
#SBATCH -t 8:00:00
#SBATCH -o ../../../processed-data/09_HD_cell_level/no_secondary/registration_banksy/logs/02_cor_heatmap_all_genes_%a.txt
#SBATCH -e ../../../processed-data/09_HD_cell_level/no_secondary/registration_banksy/logs/02_cor_heatmap_all_genes_%a.txt
#SBATCH --array=3

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

Rscript 02_cor_heatmap.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/

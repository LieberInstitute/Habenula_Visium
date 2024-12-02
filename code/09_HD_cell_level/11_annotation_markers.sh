#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=11_annotation_markers
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../processed-data/09_HD_cell_level/logs/11_annotation_markers_%a.txt
#SBATCH -e ../../processed-data/09_HD_cell_level/logs/11_annotation_markers_%a.txt
#SBATCH --array=2-28%10

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

Rscript 11_annotation_markers.R

echo "**** Job ends ****"
date

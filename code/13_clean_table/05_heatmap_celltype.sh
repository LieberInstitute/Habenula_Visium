#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=15G
#SBATCH --job-name=05_heatmap_celltype
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../processed-data/13_clean_table/logs/05_heatmap_celltype.txt
#SBATCH -e ../../processed-data/13_clean_table/logs/05_heatmap_celltype.txt

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
module load liana_plus

## List current modules for reproducibility
module list

python3 05_heatmap_celltype_new.py

echo "**** Job ends ****"
date
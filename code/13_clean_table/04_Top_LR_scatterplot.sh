#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=15G
#SBATCH --job-name=04_Top_LR_scatterplot
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../processed-data/13_clean_table/logs/04_Top_LR_scatterplot.txt
#SBATCH -e ../../processed-data/13_clean_table/logs/04_Top_LR_scatterplot.txt

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

python3 04_Top_LR_scatterplot_new.py

echo "**** Job ends ****"
date
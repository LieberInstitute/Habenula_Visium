#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=20G
#SBATCH --job-name=25_ficture_manual_plot
#SBATCH -o ../../../processed-data/10_HD_bin_level/no_secondary/cell_environment/logs/25_ficture_manual_plot_%a.txt
#SBATCH -e ../../../processed-data/10_HD_bin_level/no_secondary/cell_environment/logs/25_ficture_manual_plot_%a.txt
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH --array=17

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module list

module load conda_R/4.5
Rscript 25_ficture_manual_plot.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/

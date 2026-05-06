#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=5G
#SBATCH --job-name=19_plot_ficture_crawdad
#SBATCH -o ../../../processed-data/10_HD_bin_level/no_secondary/cell_environment/logs/19_plot_ficture_crawdad_%a.txt
#SBATCH -e ../../../processed-data/10_HD_bin_level/no_secondary/cell_environment/logs/19_plot_ficture_crawdad_%a.txt
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH --array=3-10,20%9

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
Rscript 19_plot_ficture_crawdad.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/

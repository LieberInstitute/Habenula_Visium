#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=25G
#SBATCH --job-name=20_extra_registration
#SBATCH -o ../../../processed-data/10_HD_bin_level/no_secondary/cell_environment/logs/20_extra_registration_%a.txt
#SBATCH -e ../../../processed-data/10_HD_bin_level/no_secondary/cell_environment/logs/20_extra_registration_%a.txt
#SBATCH -c 4
#SBATCH -t 1-0:00:00
#SBATCH --array=3
#SBATCH --reservation=neagles-2wk

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
Rscript 20_extra_registration.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/

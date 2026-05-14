#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=80G
#SBATCH --job-name=05_registration_inputs
#SBATCH -o ../../../../processed-data/10_HD_bin_level/no_secondary/cell_environment/investigation/logs/05_registration_inputs_%a.txt
#SBATCH -e ../../../../processed-data/10_HD_bin_level/no_secondary/cell_environment/investigation/logs/05_registration_inputs_%a.txt
#SBATCH -c 8
#SBATCH -t 1-0:00:00
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

module list

module load conda_R/4.5
Rscript 05_registration_inputs.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/

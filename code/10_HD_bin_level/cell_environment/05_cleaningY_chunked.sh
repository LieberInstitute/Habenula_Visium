#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=40G
#SBATCH --job-name=05_cleaningY_chunked
#SBATCH -o ../../../processed-data/10_HD_bin_level/probe_fix/cell_environment/logs/05_cleaningY_chunked.txt
#SBATCH -e ../../../processed-data/10_HD_bin_level/probe_fix/cell_environment/logs/05_cleaningY_chunked.txt
#SBATCH -c 1
#SBATCH -t 8:00:00
#SBATCH --array=1-50%20

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

module load conda_R/4.4.x
Rscript 05_cleaningY_chunked.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/

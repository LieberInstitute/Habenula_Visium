#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=32G
#SBATCH --job-name=01_explore_distances
#SBATCH -o ../../../processed-data/10_HD_bin_level/probe_fix/cell_environment/logs/01_explore_distances_%a.txt
#SBATCH -e ../../../processed-data/10_HD_bin_level/probe_fix/cell_environment/logs/01_explore_distances_%a.txt
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH --array=1-5%5

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

module load visium_hd/1.0
python 01_explore_distances.py

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/

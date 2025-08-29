#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=14_2_LIANA+_preprocess
#SBATCH -c 1
#SBATCH -t 2:00:00
#SBATCH -o ../../../processed-data/10_HD_bin_level/logs/14_2_LIANA+_preprocess_%a.txt
#SBATCH -e ../../../processed-data/10_HD_bin_level/logs/14_2_LIANA+_preprocess_%a.txt
#SBATCH --array=1-3%3

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
module load visium_hd/1.0

## List current modules for reproducibility
module list

python3 14_2_LIANA+_preprocess.py

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/

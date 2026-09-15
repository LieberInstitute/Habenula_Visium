#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=100G
#SBATCH --job-name=1_preprocess
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../../../processed-data/10_HD_bin_level/no_secondary/liana2/logs/1_preprocess.txt
#SBATCH -e ../../../../processed-data/10_HD_bin_level/no_secondary/liana2/logs/1_preprocess.txt

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

python 1_preprocess.py
python 2_LIANA+_preprocess_cell.py

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/

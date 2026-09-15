#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=02_Disease_LRs
#SBATCH -t 1-00:00:00
#SBATCH -o ../../../../processed-data/10_HD_bin_level/no_secondary/liana2/open_targets/logs/02_Disease_LRs.txt
#SBATCH -e ../../../../processed-data/10_HD_bin_level/no_secondary/liana2/open_targets/logs/02_Disease_LRs.txt

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load the python module
module load liana_plus/1.7.1

## List current modules for reproducibility
module list

## Edit with your job command
python 02_Disease_LRs.py --disease Substance_dependence
python 02_Disease_LRs.py --disease MDD

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.4.0
## available from http://research.libd.org/slurmjobs/

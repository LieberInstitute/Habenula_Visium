#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=100G
#SBATCH --job-name=19_1_LIANA+_cell_bandwidth
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../../processed-data/10_HD_bin_level/no_secondary/liana2/logs/19_1_LIANA+_cell_bandwidth_%a.txt
#SBATCH -e ../../../processed-data/10_HD_bin_level/no_secondary/liana2/logs/19_1_LIANA+_cell_bandwidth_%a.txt
#SBATCH --array=1-8%8

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
module load liana_plus/1.7.1

## List current modules for reproducibility
module list

bandwidth=7500

python3 '19_1_LIANA+_cell_bandwidth.py' --bandwidth "$bandwidth"

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/

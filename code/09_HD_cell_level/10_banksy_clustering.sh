#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=10_banksy_clustering
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --array=1-20%20

## Define loops and appropriately subset each variable for the array task ID
all_res=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1)
res=${all_res[$(( $SLURM_ARRAY_TASK_ID / 2 % 10 ))]}

all_lambda=(0.2 0.8)
lambda=${all_lambda[$(( $SLURM_ARRAY_TASK_ID / 1 % 2 ))]}

## Explicitly pipe script output to a log
log_path=../../processed-data/09_HD_cell_level/logs/10_banksy_clustering_${res}_${lambda}_${SLURM_ARRAY_TASK_ID}.txt

{
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
module load conda_R/4.4.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 10_banksy_clustering.R --res ${res} --lambda ${lambda}

echo "**** Job ends ****"
date

} > $log_path 2>&1

## This script was made using slurmjobs version 1.3.0
## available from http://research.libd.org/slurmjobs/

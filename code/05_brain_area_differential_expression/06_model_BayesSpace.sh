#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=30G
#SBATCH --job-name=06_model_BayesSpace
#SBATCH -c 2
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
# SBATCH --mail-type=ALL
#SBATCH --array=2-28%20


## Define loops and appropriately subset each variable for the array task ID

#all_selected_BS_k=(13 21 26)
# all_selected_BS_k=($(seq 2 28))
# BS_k=${all_selected_BS_k[$(( $SLURM_ARRAY_TASK_ID / 1 % 28 ))]}
# 13 = 1 Hb domain
# 21 = 2 Hb domains
# 26 = 3 Hb domains

## Explicitly pipe script output to a log
log_path=logs/06_model_BayesSpace_${SLURM_ARRAY_TASK_ID}.txt

{
set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
echo "SLURM_ARRAY_TASK_ID is: $SLURM_ARRAY_TASK_ID"

## Load the R module
module load conda_R/4.4.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 06_model_BayesSpace.R
#Rscript 06_model_BayesSpace.R --BS_k ${BS_k}

echo "**** Job ends ****"
date


} > $log_path 2>&1

## This script was made using slurmjobs version 1.2.1
## available from http://research.libd.org/slurmjobs/

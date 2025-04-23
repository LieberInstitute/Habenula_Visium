#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=30G
#SBATCH --job-name=02_explore_expr_variability
#SBATCH -c 2
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
# SBATCH --mail-type=ALL
#SBATCH --array=1-4%4

## Define BayesSpace k of interest
#BS_k=(3 9 17)
BS_k=(3 13 21 26)
BS_k=${BS_k[$(( $SLURM_ARRAY_TASK_ID / 1 % 4 ))]}

## Explicitly pipe script output to a log
log_path=logs/02_explore_expr_variability_BS${BS_k}.txt

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

## Load the R module
module load conda_R/4.4.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 02_explore_expr_variability.R --BS_k ${BS_k}

echo "**** Job ends ****"
date

  
} > $log_path 2>&1

## This script was made using slurmjobs version 1.2.1
## available from http://research.libd.org/slurmjobs/

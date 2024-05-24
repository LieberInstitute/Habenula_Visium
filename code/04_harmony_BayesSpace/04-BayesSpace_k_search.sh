#!/bin/bash
#SBATCH -p bluejay
#SBATCH --mem=32G
#SBATCH --job-name=BayesSpace_k_search
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o logs/03-BayesSpace_k_search.%a.txt
#SBATCH -e logs/03-BayesSpace_k_search.%a.txt
#SBATCH --mail-type=ALL
#SBATCH --array=2-28%20

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
module load conda_R/4.3.x

## List current modules for reproducibility
module list

Rscript 03-BayesSpace_k_search.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.0
## available from http://research.libd.org/slurmjobs/

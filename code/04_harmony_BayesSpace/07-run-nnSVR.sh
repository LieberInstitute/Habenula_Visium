#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=20G
#SBATCH --job-name=05-run_nnSVG
#SBATCH -o ../../processed-data/05_harmony_BayesSpace/logs/05-run_nnSVG_%a.log
#SBATCH -e ../../processed-data/05_harmony_BayesSpace/logs/05-run_nnSVG_%a.log
#SBATCH --array=1-10%5

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load conda_R/4.3
Rscript 05-run_nnSVG.R

echo "**** Job ends ****"
date
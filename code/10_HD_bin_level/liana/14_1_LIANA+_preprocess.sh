#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=32G
#SBATCH --job-name=14_1_LIANA+_preprocess
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../../processed-data/10_HD_bin_level/new_samples2/liana/logs/14_1_LIANA+_preprocess.txt
#SBATCH -e ../../../processed-data/10_HD_bin_level/new_samples2/liana/logs/14_1_LIANA+_preprocess.txt

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load visium_hd/1.0
module list
module load conda_R/4.4

Rscript 14_1_LIANA+_preprocess.R

echo "**** Job ends ****"
date

#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=100G
#SBATCH --job-name=08_batcheffect
#SBATCH -c 1
#SBATCH -t 96:00:00
#SBATCH -o ../../../processed-data/10_HD_bin_level/ficture_harmony/logs/08_batcheffect.txt
#SBATCH -e ../../../processed-data/10_HD_bin_level/ficture_harmony/logs/08_batcheffect.txt

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
module load conda_R

## List current modules for reproducibility
module list

Rscript 08_batcheffect.R

echo "**** Job ends ****"
date

#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=80G
#SBATCH --job-name=12_preprocess_anndata
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../../processed-data/10_HD_bin_level/liana/logs/12_preprocess_anndata.txt
#SBATCH -e ../../../processed-data/10_HD_bin_level/liana/logs/12_preprocess_anndata.txt
#SBATCH --array=1-5%5

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

## List current modules for reproducibility
module list

python 12_preprocess_anndata.py

echo "**** Job ends ****"
date

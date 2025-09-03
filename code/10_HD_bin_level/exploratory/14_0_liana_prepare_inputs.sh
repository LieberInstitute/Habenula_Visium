#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=32G
#SBATCH --job-name=14_0_liana_prepare_inputs
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../../processed-data/09_HD_cell_level/probe_fix/logs/14_0_liana_prepare_inputs.txt
#SBATCH -e ../../../processed-data/09_HD_cell_level/probe_fix/logs/14_0_liana_prepare_inputs.txt

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

python 14_0_liana_prepare_inputs.py

echo "**** Job ends ****"
date

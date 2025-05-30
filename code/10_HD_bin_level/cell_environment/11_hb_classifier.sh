#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=11_hb_classifier
#SBATCH -o ../../../processed-data/10_HD_bin_level/probe_fix/cell_environment/logs/11_hb_classifier_%a.txt
#SBATCH -e ../../../processed-data/10_HD_bin_level/probe_fix/cell_environment/logs/11_hb_classifier_%a.txt
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH --array=20

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module list

module load visium_hd/1.0
python 11_hb_classifier.py

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/

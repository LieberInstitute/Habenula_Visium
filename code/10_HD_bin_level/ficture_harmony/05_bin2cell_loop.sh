#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=80G
#SBATCH --job-name=05_bin2cell
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../../processed-data/10_HD_bin_level/ficture_harmony/logs/05_bin2cell_loop_%a.txt
#SBATCH -e ../../../processed-data/10_HD_bin_level/ficture_harmony/logs/05_bin2cell_loop_%a.txt
#SBATCH --array=2-25%24

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

MY_VAR=${SLURM_ARRAY_TASK_ID} sbatch --export=ALL,MY_VAR 05_bin2cell_new.sh

echo "**** Job ends ****"
date
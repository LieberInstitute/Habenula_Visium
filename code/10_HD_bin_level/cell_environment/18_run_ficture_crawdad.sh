#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=5G
#SBATCH --job-name=18_run_ficture_crawdad
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH --array=1-72%20

## Define loops and appropriately subset each variable for the array task ID
all_sample_id=(Br9090_1 Br9090_2 Br8433_1 Br8433_2 Br8667_1 Br8667_2 Br3942_1 Br3942_2)
sample_id=${all_sample_id[$(( $SLURM_ARRAY_TASK_ID / 9 % 8 ))]}

all_k=(3 4 5 6 7 8 9 10 20)
k=${all_k[$(( $SLURM_ARRAY_TASK_ID / 1 % 9 ))]}

## Explicitly pipe script output to a log
log_path=../../../processed-data/10_HD_bin_level/no_secondary/cell_environment/logs/18_run_ficture_crawdad_${sample_id}_${k}_${SLURM_ARRAY_TASK_ID}.txt

{
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
module load conda_R/4.5

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 18_run_ficture_crawdad.R --sample_id ${sample_id} --k ${k}

echo "**** Job ends ****"
date

} > $log_path 2>&1

## This script was made using slurmjobs version 1.3.0
## available from http://research.libd.org/slurmjobs/



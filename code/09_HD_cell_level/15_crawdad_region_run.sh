#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=5G
#SBATCH --job-name=15_crawdad_region_run
#SBATCH -c 4
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --array=9

## Define loops and appropriately subset each variable for the array task ID
all_sample_id=(H1-W369TJK_D1_9090 H1-MVPY9BW_A1_8433 H1-MVPY9BW_D1_8667 H1-6FX4YN3_A1_3942 H1-6FX4YN3_D1_9902)
sample_id=${all_sample_id[$(( $SLURM_ARRAY_TASK_ID / 2 % 5 ))]}

all_region=(habenula thalamus)
region=${all_region[$(( $SLURM_ARRAY_TASK_ID / 1 % 2 ))]}

## Explicitly pipe script output to a log
log_path=../../processed-data/09_HD_cell_level/new_samples2/logs/15_crawdad_region_run_${sample_id}_${region}_${SLURM_ARRAY_TASK_ID}.txt

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
Rscript 15_crawdad_region_run.R --sample_id ${sample_id} --region ${region}

echo "**** Job ends ****"
date

} > $log_path 2>&1

## This script was made using slurmjobs version 1.3.0
## available from http://research.libd.org/slurmjobs/


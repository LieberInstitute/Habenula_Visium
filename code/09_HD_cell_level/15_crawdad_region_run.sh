#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=5G
#SBATCH --job-name=15_crawdad_region_run
#SBATCH -c 4
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --array=1-20%10

## Define loops and appropriately subset each variable for the array task ID
all_sample_id=(Br9090_1 Br9090_2 Br8433_1 Br8433_2 Br8667_1 Br8667_2 Br3942_1 Br3942_2 Br9902_1 Br9902_2)
sample_id=${all_sample_id[$(( $SLURM_ARRAY_TASK_ID / 2 % 10 ))]}

all_region=(habenula thalamus)
region=${all_region[$(( $SLURM_ARRAY_TASK_ID / 1 % 2 ))]}

## Explicitly pipe script output to a log
log_path=../../processed-data/09_HD_cell_level/no_secondary/logs/15_crawdad_region_run_${sample_id}_${region}_${SLURM_ARRAY_TASK_ID}.txt

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


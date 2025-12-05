#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=5G
#SBATCH --job-name=14_deciding_k
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --array=1-4%4

## Define loops and appropriately subset each variable for the array task ID
all_ref_name=(snRNAseq_fine snRNAseq_broad)
ref_name=${all_ref_name[$(( $SLURM_ARRAY_TASK_ID / 2 % 2 ))]}

all_input_method=(normalized cleaning_y)
input_method=${all_input_method[$(( $SLURM_ARRAY_TASK_ID / 1 % 2 ))]}

## Explicitly pipe script output to a log
log_path=../../../processed-data/10_HD_bin_level/new_samples2/ficture_harmony/logs/14_deciding_k_${ref_name}_${input_method}_${SLURM_ARRAY_TASK_ID}.txt

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
module load conda_R/4.4

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 14_deciding_k.R --ref_name ${ref_name} --input_method ${input_method}

echo "**** Job ends ****"
date

} > $log_path 2>&1

## This script was made using slurmjobs version 1.3.0
## available from http://research.libd.org/slurmjobs/

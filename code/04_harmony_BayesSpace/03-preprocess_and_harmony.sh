#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=100G
#SBATCH --job-name=03-preprocess_and_harmony
#SBATCH -c 2
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
# SBATCH --mail-type=ALL

## Explicitly pipe script output to a log
log_path=logs/03-preprocess_and_harmony.txt

# if [[ ! -z $SLURMD_NODENAME ]]; then
#     job_id=$SLURM_JOB_ID
#     job_name=$SLURM_JOB_NAME
#     node_name=$SLURMD_NODENAME
# else
#     job_id=$JOB_ID
#     job_name=$JOB_NAME
#     node_name=$HOSTNAME
# fi

{
set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${job_id}"
echo "Job name: ${job_name}"
echo "Node name: ${node_name}"

## List current modules for reproducibility
module load conda_R/4.4.x
module list

Rscript 03-preprocess_and_harmony.R

echo "**** Job ends ****"
date

} > $log_path 2>&1

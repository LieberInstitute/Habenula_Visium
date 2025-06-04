#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=15G
#SBATCH --job-name=09_comparing_BayesSpace_clusters
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
# SBATCH --mail-type=ALL


## Explicitly pipe script output to a log
# mkdir -p logs
log_path=logs/09_comparing_BayesSpace_clusters_${SLURM_ARRAY_TASK_ID}.txt

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
module load conda_R/4.4.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 09_comparing_BayesSpace_clusters.R

echo "**** Job ends ****"
date
echo "Exit code: $ret"
exit $ret

} > $log_path 2>&1

## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/

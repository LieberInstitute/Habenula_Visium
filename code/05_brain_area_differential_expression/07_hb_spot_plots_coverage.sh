#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=30G
#SBATCH --job-name=07_hb_spot_plots_coverage
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --array=1-4%4
# SBATCH --mail-type=ALL

## Define BayesSpace k of interest
BS_k_list=(3 11 15 20)
BS_k=${BS_k_list[$((SLURM_ARRAY_TASK_ID - 1))]}

mkdir -p logs
log_path=logs/07_hb_spot_plots_coverage${BS_k}.txt

{
set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
echo "BS_k value: ${BS_k}"

## Load the R module
module load conda_R/4.4.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 07_hb_spot_plots_coverage.R --BS_k ${BS_k}
ret=$?

echo "**** Job ends ****"
date
echo "Exit code: $ret"
exit $ret

} > $log_path 2>&1

## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/
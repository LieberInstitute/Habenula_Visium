#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=2G
#SBATCH --job-name=05_bin2cell_wrapper
#SBATCH -c 1
#SBATCH -t 10:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --array=2-25%2

export OUTER_TASK=$SLURM_ARRAY_TASK_ID
log_path=../../../processed-data/10_HD_bin_level/ficture_harmony/logs/05_bin2cell_k${OUTER_TASK}_sample%a.txt
sbatch \
    --export=ALL,OUTER_TASK \
    -o $log_path \
    -e $log_path \
    05_bin2cell.sh

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/

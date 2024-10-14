#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=32G
#SBATCH --job-name=01_build_spe
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../processed-data/10_HD_bin_level/logs/01_build_spe.txt
#SBATCH -e ../../processed-data/10_HD_bin_level/logs/01_build_spe.txt

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

repo_dir=$(git rev-parse --show-toplevel)
spatial_dir=$repo_dir/processed-data/01_spaceranger/H1-W369TJK_D1_9090/outs/binned_outputs/square_008um/spatial

#   Get spatial coordinates as a CSV
if [[ ! -f $spatial_dir/tissue_positions.parquet ]]; then
    echo "Converting spatial coords to CSV..."
    module load ficture/0.0.3.1
    parquet-tools csv $spatial_dir/tissue_positions.parquet \
        > $spatial_dir/tissue_positions.csv
    module unload ficture
fi

## Load the R module
module load conda_R/4.4

## List current modules for reproducibility
module list

Rscript 01_build_spe.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/

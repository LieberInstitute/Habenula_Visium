#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=64G
#SBATCH --job-name=01_build_spe
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../../processed-data/10_HD_bin_level/new_samples2/ficture_harmony/logs/01_build_spe.txt
#SBATCH -e ../../../processed-data/10_HD_bin_level/new_samples2/ficture_harmony/logs/01_build_spe.txt

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
sample_info_path=$repo_dir/raw-data/sample_info/hd_basic_info.csv

#   Get spatial coordinates as a CSV (from parquet format) for each sample where
#   it doesn't exist
module load visium_hd/1.0
for i in $(seq 1 5); do
    sample_id=$(awk -F, "NR==${i}+1 {print \$1}" $sample_info_path)
    spatial_dir=$repo_dir/$(awk -F, "NR==${i}+1 {print \$2}" $sample_info_path)/outs/binned_outputs/square_002um/spatial

    if [[ ! -f $spatial_dir/tissue_positions.csv ]]; then
        echo "Converting spatial coords to CSV for sample ${sample_id}..."
        
        parquet-tools csv $spatial_dir/tissue_positions.parquet \
            > $spatial_dir/tissue_positions.csv
    fi
done
module unload visium_hd

## Load the R module
module load conda_R/4.5

## List current modules for reproducibility
module list

Rscript 01_build_spe.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/

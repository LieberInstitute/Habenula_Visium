#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=200G
#SBATCH --job-name=07_ficture_run
#SBATCH -c 1
#SBATCH -t 5-0:00:00
#SBATCH -o ../../../processed-data/10_HD_bin_level/no_secondary/cell_environment/logs/07_ficture_run_%a.txt
#SBATCH -e ../../../processed-data/10_HD_bin_level/no_secondary/cell_environment/logs/07_ficture_run_%a.txt
#SBATCH --array=3,5%2

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load visium_hd/1.0

#   Path definitions
repo_dir=$(git rev-parse --show-toplevel)
in_dir=$repo_dir/processed-data/10_HD_bin_level/no_secondary/cell_environment/ficture_inputs
out_dir=$repo_dir/processed-data/10_HD_bin_level/no_secondary/cell_environment/ficture_outputs/cleaningy/k_${SLURM_ARRAY_TASK_ID}

mkdir -p $out_dir

#   Run full FICTURE pipeline
ficture run_together \
    --in-tsv $in_dir/cleaningy_input.tsv.gz \
    --in-minmax $in_dir/cleaningy_minmax.tsv \
    --out-dir $out_dir \
    --mu-scale 1 \
    --major-axis X \
    --all \
    --n-factor ${SLURM_ARRAY_TASK_ID} \
    --fractional-count 1

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.2
## available from http://research.libd.org/slurmjobs/

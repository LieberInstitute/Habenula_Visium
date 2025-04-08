#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=100G
#SBATCH --job-name=03_ficture_run
#SBATCH -c 1
#SBATCH -t 3-0:00:00
#SBATCH -o ../../../processed-data/10_HD_bin_level/logs/all_sample_002um/03_ficture_run_%a.log
#SBATCH -e ../../../processed-data/10_HD_bin_level/logs/all_sample_002um/03_ficture_run_%a.log
#SBATCH --array=3,7-9,20-25%5

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

ml ficture/dev_a455e5c

repo_dir=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium

#   Path definitions
in_dir=$repo_dir/processed-data/10_HD_bin_level/ficture/inputs/all_samples
out_dir=$repo_dir/processed-data/10_HD_bin_level/ficture/outputs/all_samples/k_$SLURM_ARRAY_TASK_ID

mkdir -p $out_dir

#   Run full FICTURE pipeline
ficture run_together \
    --in-tsv $in_dir/transcripts.merged.sorted.tsv.gz \
    --in-minmax $in_dir/merged_minmax.tsv \
    --out-dir $out_dir \
    --mu-scale 1 \
    --major-axis X \
    --all \
    --n-factor ${SLURM_ARRAY_TASK_ID}

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.2
## available from http://research.libd.org/slurmjobs/

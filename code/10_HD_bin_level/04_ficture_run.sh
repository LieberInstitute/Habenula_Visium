#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=50G
#SBATCH --job-name=04_ficture_run
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../processed-data/10_HD_bin_level/logs/04_ficture_run_%a.log
#SBATCH -e ../../processed-data/10_HD_bin_level/logs/04_ficture_run_%a.log
#SBATCH --array=1-5%5

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

repo_dir=$(git rev-parse --show-toplevel)
sample_id_path=$repo_dir/raw-data/sample_info/hd_sample_list.txt

this_sample=$(awk "NR==${SLURM_ARRAY_TASK_ID}" $sample_id_path)

#   Path definitions
in_dir=$repo_dir/processed-data/10_HD_bin_level/ficture/inputs/$this_sample
out_dir=$repo_dir/processed-data/10_HD_bin_level/ficture/outputs/$this_sample

mkdir -p $out_dir

#   Run full FICTURE pipeline
ficture run_together \
    --in-tsv $in_dir/transcripts.sorted.tsv.gz \
    --in-minmax $in_dir/minmax.tsv \
    --out-dir $out_dir \
    --mu-scale 1 \
    --major-axis X \
    --all

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.2
## available from http://research.libd.org/slurmjobs/

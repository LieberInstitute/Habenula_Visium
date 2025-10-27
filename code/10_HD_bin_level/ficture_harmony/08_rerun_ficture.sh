#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=200G
#SBATCH --job-name=08_rerun_ficture
#SBATCH -c 1
#SBATCH -t 2-0:00:00
#SBATCH -o ../../../processed-data/10_HD_bin_level/new_samples/ficture_harmony/logs/08_rerun_ficture_%a.txt
#SBATCH -e ../../../processed-data/10_HD_bin_level/new_samples/ficture_harmony/logs/08_rerun_ficture_%a.txt
#SBATCH --array=3-40,70,100%10

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
in_dir=$repo_dir/processed-data/10_HD_bin_level/new_samples/ficture_harmony/ficture_inputs/cleany
out_dir=$repo_dir/processed-data/10_HD_bin_level/new_samples/ficture_harmony/ficture_outputs/cleany/k_${SLURM_ARRAY_TASK_ID}

mkdir -p $out_dir

#   Run full FICTURE pipeline
if [[ $SLURM_ARRAY_TASK_ID -eq 2 ]]; then
    #   Overwride the default of finding top 3 factors, since only 2 exist
    ficture run_together \
        --in-tsv $in_dir/input.tsv.gz \
        --in-minmax $in_dir/minmax.tsv \
        --out-dir $out_dir \
        --mu-scale 1 \
        --major-axis X \
        --all \
        --n-factor ${SLURM_ARRAY_TASK_ID} \
        --decode-top-k 2 \
        --fractional-count 1
else
    ficture run_together \
        --in-tsv $in_dir/input.tsv.gz \
        --in-minmax $in_dir/minmax.tsv \
        --out-dir $out_dir \
        --mu-scale 1 \
        --major-axis X \
        --all \
        --n-factor ${SLURM_ARRAY_TASK_ID} \
        --fractional-count 1
fi

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.2
## available from http://research.libd.org/slurmjobs/

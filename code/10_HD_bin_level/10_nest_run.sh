#!/bin/bash
#SBATCH -p caracol
#SBATCH --mem=32G
#SBATCH --job-name=10_nest_run
#SBATCH -c 1
#SBATCH -t 5-0:00:00
#SBATCH -o ../../processed-data/10_HD_bin_level/logs/10_nest_run.txt
#SBATCH -e ../../processed-data/10_HD_bin_level/logs/10_nest_run.txt

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

## List current modules for reproducibility
module list

this_sample=H1-W369TJK_D1_9090
repo_dir=$(git rev-parse --show-toplevel)
in_dir=$repo_dir/processed-data/01_spaceranger/$this_sample/outs/binned_outputs/square_008um
out_dir=$repo_dir/processed-data/10_HD_bin_level/nest
num_runs=5

#   While ordinarily this might be parallelized as an array job, in practice
#   GPU availability is so scarce that it's faster to hold on to a GPU and
#   run serially
for i in $(seq 1 $num_runs); do
    echo "---- Starting iteration $i ----"
    nest run \
        --data_name=$this_sample \
        --metadata_to=$out_dir \
        --model_path=$out_dir \
        --embedding_path=$out_dir \
        --manual_seed=yes \
        --seed=$i \
        --model_name=NEST_model \
        --run_id=$i
done

echo "**** Job ends ****"
date

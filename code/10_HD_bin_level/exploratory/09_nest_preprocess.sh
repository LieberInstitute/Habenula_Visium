#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=200G
#SBATCH --job-name=09_nest_preprocess
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH -o ../../processed-data/10_HD_bin_level/logs/09_nest_preprocess.txt
#SBATCH -e ../../processed-data/10_HD_bin_level/logs/09_nest_preprocess.txt

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
db_path=/jhpce/shared/libd/core/visium_hd/1.0/NEST/database/NEST_database.csv

mkdir -p $out_dir

nest preprocess \
    --data_name=$this_sample \
    --data_from=$in_dir \
    --database_path=$db_path \
    --data_to=$out_dir \
    --metadata_to=$out_dir

echo "**** Job ends ****"
date

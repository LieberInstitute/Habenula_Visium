#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=2G
#SBATCH --job-name=02_link_shiny_files
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../processed-data/12_apps_and_sharing/logs/02_link_shiny_files.txt
#SBATCH -e ../../processed-data/12_apps_and_sharing/logs/02_link_shiny_files.txt

# Symlink Shiny app files so relative paths can be used from the destination
# directory for the app

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
source_dir=${repo_dir}/processed-data/12_apps_and_sharing/01_prep_objects
dest_dir=${repo_dir}/processed-data/12_apps_and_sharing/shiny_app

modeling_path=${repo_dir}/processed-data/09_HD_cell_level/no_secondary/registration_banksy/modeling_results/1_8_cell_types.rds

for f_name in spe_shiny.qs2 spe_pb_shiny.qs2 sig_genes_shiny.qs2; do
    rm -f ${dest_dir}/${f_name}
    ln -s ${source_dir}/${f_name} ${dest_dir}/${f_name}
done

rm -f ${dest_dir}/modeling_results.rds
ln -s ${modeling_path} ${dest_dir}/modeling_results.rds

echo "**** Job ends ****"
date

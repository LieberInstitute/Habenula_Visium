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
dest_dir=${repo_dir}/code/12_apps_and_sharing/shiny_app

#   Using relative paths for sym links so things work with git
cd ${dest_dir}
for f_name in spe_shiny.qs2 sce_pb_shiny.qs2 sig_genes_shiny.qs2; do
    rm -f ${f_name}
    ln -s ../../../processed-data/12_apps_and_sharing/01_prep_objects/${f_name} ${f_name}
done

rm -f modeling_results.rds
ln -s ../../../processed-data/09_HD_cell_level/no_secondary/registration_banksy/modeling_results/1_8_cell_types.rds modeling_results.rds

echo "**** Job ends ****"
date

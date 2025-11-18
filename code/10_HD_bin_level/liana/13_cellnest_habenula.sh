#!/bin/bash
#SBATCH -p caracol
#SBATCH --mem=32G
#SBATCH --job-name=13_cellnest_habenula
#SBATCH -c 1
#SBATCH --gres=gpu:1
#SBATCH -t 1:00:00
#SBATCH -o ../../../processed-data/10_HD_bin_level/logs/13_cellnest_habenula.txt
#SBATCH -e ../../../processed-data/10_HD_bin_level/logs/13_cellnest_habenula.txt

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
in_dir=$repo_dir/processed-data/10_HD_bin_level/nest/input_habenula/${this_sample}.h5ad
out_dir=$repo_dir/processed-data/10_HD_bin_level/nest

module load singularity
# navigate to the directory that has CellNEST repository
cd /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/nest/CellNEST/
echo "Current working directory: `pwd`"

singularity shell --home=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/nest/ /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/nest/cellnest_container/cellnest_image.sif
bash cellnest preprocess --data_name="${this_sample}_node_spatial" --data_type=anndata --data_from="${in_dir}" --split=1

exit

echo "**** Job ends ****"
date
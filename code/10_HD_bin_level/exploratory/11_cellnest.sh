#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=500G
#SBATCH --job-name=11_cellnest_preprocess
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH -o ../../../processed-data/10_HD_bin_level/logs/11_cellnest_preprocess.txt
#SBATCH -e ../../../processed-data/10_HD_bin_level/logs/11_cellnest_preprocess.txt

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
in_dir=$repo_dir/processed-data/01_spaceranger/$this_sample/outs/binned_outputs/square_002um
out_dir=$repo_dir/processed-data/10_HD_bin_level/nest


module load singularity
# navigate to the directory that has CellNEST repository
cd /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/nest/CellNEST/
echo "Current working directory: `pwd`"

# python to_anndata.py \
#   --data_from=$in_dir \
#   --data_to_path=$in_dir \
#   --file_name=$this_sample \
#   --tissue_position_file=data/LUAD_GSM5702473_TD1/GSM5702473_TD1_tissue_positions_list.csv

singularity shell --home=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/nest/ /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/nest/cellnest_container/cellnest_image.sif
bash cellnest preprocess --data_name="${this_sample}_node_spatial" --data_type=visium --data_from="${in_dir}" --split=1 --filter_min_cell=3
# bash cellnest_visium preprocess --data_name="${this_sample}_node_spatial" --data_from="${in_dir}" --split=1

exit

echo "**** Job ends ****"
date
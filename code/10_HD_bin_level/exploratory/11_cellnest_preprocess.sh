# Preprocess the CellNEST data for HD bin level analysis
# One sample



#singularity pull cellnest_image.sif library://fatema/collection/cellnest_image.sif:latest
singularity shell  --nv --home=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/nest/ /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/nest/cellnest_container/cellnest_image.sif
Singularity> bash cellnest preprocess --data_name='V1_Human_Lymph_Node_spatial' --data_from='/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/nest/CellNEST/data/V1_Human_Lymph_Node_spatial/'

nohup bash cellnest run --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --model_name='CellNEST_V1_Human_Lymph_Node_spatial' --run_id=1 > output_human_lymph_node_run1.log &


Singularity> bash cellnest preprocess --data_name='V1_Human_Lymph_Node_spatial' --data_from='/cluster/projects/prof-group/fatema/CellNEST/data/V1_Human_Lymph_Node_spatial/'
Singularity> ls metadata/V1_Human_Lymph_Node_spatial/
Singularity> ls input_graph/V1_Human_Lymph_Node_spatial/
$ singularity shell --nv --home=/cluster/projects/prof-group/fatema/cellnest_container/ /cluster/projects/prof-group/fatema/cellnest_container/cellnest_image.sif


repo_dir=$(git rev-parse --show-toplevel)
out_dir=$repo_dir/processed-data/10_HD_bin_level/nest
mkdir cellnest_container 
cd cellnest_container


#!/bin/bash
# ---------------------------------------------------------------------
# SLURM script for a job on a cluster.
# ---------------------------------------------------------------------
#SBATCH -c 16
#SBATCH --mem=30GB
#SBATCH --time=72:00:00
#SBATCH --job-name=CellNEST_lymph
#SBATCH --output=some_name-%j.out
# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
# ---------------------------------------------------------------------
echo ""
echo "Job Array ID / Job ID: $SLURM_ARRAY_JOB_ID / $SLURM_JOB_ID"
echo "This is job $SLURM_ARRAY_TASK_ID out of $SLURM_ARRAY_TASK_COUNT jobs."
echo ""
# ---------------------------------------------------------------------

# module load
module load singularity

# navigate to the directory that has CellNEST repository
cd /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/nest/CellNEST/
echo "Current working directory: `pwd`"

# cellnest_image.sif is here: /cluster/projects/prof-group/fatema/cellnest_container/cellnest_image.sif

# Preprocess the CellNEST data
singularity shell  --nv --home=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/nest/ /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/nest/cellnest_container/cellnest_image.sif
bash cellnest preprocess --data_name='V1_Human_Lymph_Node_spatial' --data_from='/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/nest/CellNEST/data/V1_Human_Lymph_Node_spatial/'

# checking the available gpu memory before starting the training
nvidia-smi

# run your python script with parameters using singularity
echo "lymph. Going to start process: run 1"

singularity run --nv --home=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/nest/ /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/nest/cellnest_container/cellnest_image.sif \
bash cellnest run --data_name='V1_Human_Lymph_Node_spatial' --data_from='/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/nest/CellNEST/data/V1_Human_Lymph_Node_spatial/' \
--num_epoch 80000 --manual_seed='yes' --seed=1 --model_name='CellNEST_V1_Human_Lymph_Node_spatial' --run_id=1

echo "lymph. Going to start process: run 2"

singularity run --nv --home=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/nest/ /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/nest/cellnest_container/cellnest_image.sif \
bash cellnest run --data_name='V1_Human_Lymph_Node_spatial' --data_from='/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/nest/CellNEST/data/V1_Human_Lymph_Node_spatial/' \
--num_epoch 80000 --manual_seed='yes' --seed=1 --model_name='CellNEST_V1_Human_Lymph_Node_spatial' --run_id=2

echo "lymph. Going to start process: run 3"

singularity run --nv --home=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/nest/ /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/nest/cellnest_container/cellnest_image.sif \
bash cellnest run --data_name='V1_Human_Lymph_Node_spatial' --data_from='/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/nest/CellNEST/data/V1_Human_Lymph_Node_spatial/' \
--num_epoch 80000 --manual_seed='yes' --seed=1 --model_name='CellNEST_V1_Human_Lymph_Node_spatial' --run_id=3

echo "lymph. Going to start process: run 4"

singularity run --nv --home=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/nest/ /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/nest/cellnest_container/cellnest_image.sif \
bash cellnest run --data_name='V1_Human_Lymph_Node_spatial' --data_from='/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/nest/CellNEST/data/V1_Human_Lymph_Node_spatial/' \
--num_epoch 80000 --manual_seed='yes' --seed=1 --model_name='CellNEST_V1_Human_Lymph_Node_spatial' --run_id=4

echo "lymph. Going to start process: run 5"

singularity run --nv --home=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/nest/ /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/nest/cellnest_container/cellnest_image.sif \
bash cellnest run --data_name='V1_Human_Lymph_Node_spatial' --data_from='/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/nest/CellNEST/data/V1_Human_Lymph_Node_spatial/' \
--num_epoch 80000 --manual_seed='yes' --seed=1 --model_name='CellNEST_V1_Human_Lymph_Node_spatial' --run_id=5

echo "Job finished with exit code $? at: `date`"

cellnest postprocess --data_name='V1_Human_Lymph_Node_spatial' --model_name='CellNEST_V1_Human_Lymph_Node_spatial' --total_runs=1

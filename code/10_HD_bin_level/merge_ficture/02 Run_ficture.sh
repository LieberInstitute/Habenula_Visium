#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=50G
#SBATCH --job-name=04_ficture_02_Run_ficture.sh
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../processed-data/10_HD_bin_level/logs/04_ficture_02_Run_ficture_%a.log
#SBATCH -e ../../processed-data/10_HD_bin_level/logs/04_ficture_02_Run_ficture_%a.log
#SBATCH --array=1-5%5

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

ml ficture/dev_a455e5c
ml spatula/f0e9936

repo_dir=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium
sample_id_path=$repo_dir/raw-data/sample_info/hd_sample_list.txt
temp_dir=$MYSCRATCH
this_sample=$(awk "NR==${SLURM_ARRAY_TASK_ID}" $sample_id_path)

#   Path definitions
in_dir=$repo_dir/processed-data/10_HD_bin_level/ficture/inputs/square_008um/$this_sample
out_dir=$repo_dir/processed-data/10_HD_bin_level/ficture/outputs/square_008um/$this_sample

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
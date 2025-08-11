#!/bin/bash
#SBATCH --mem=80G
#SBATCH -c 8
#SBATCH -p katun
#SBATCH --job-name=05-hd_pilot_segment
#SBATCH -o logs/05-hd_pilot_segment_%a.txt
#SBATCH -e logs/05-hd_pilot_segment_%a.txt
#SBATCH --array=6-9%2

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## load SpaceRanger
module load spaceranger/4.0.1

## List current modules for reproducibility
module list

repo_dir=$(git rev-parse --show-toplevel)
SAMPLE=$(awk "NR==${SLURM_ARRAY_TASK_ID}" pilot_samples.txt)
IMG_PATH=$repo_dir/raw-data/images/vis-hd/pilot/${SAMPLE}.tif
OUT_DIR=$repo_dir/processed-data/01_spaceranger/probe_fix/segmentation/pilot/${SAMPLE}

mkdir -p ${OUT_DIR}

echo "Processing sample ${SAMPLE}"

spaceranger segment \
    --id=${SAMPLE} \
    --tissue-image=${IMG_PATH} \
    --output-dir=${OUT_DIR} \
    --localcores=8 \
    --localmem=80 \
    --disable-ui

echo "**** Job ends ****"
date

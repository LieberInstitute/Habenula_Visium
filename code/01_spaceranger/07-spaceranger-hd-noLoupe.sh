#!/bin/bash
#SBATCH --mem=80G
#SBATCH -c 8
#SBATCH -p katun
#SBATCH --job-name=07-spaceranger-hd-noLoupe
#SBATCH -o logs/07-spaceranger-hd-noLoupe_%a.txt
#SBATCH -e logs/07-spaceranger-hd-noLoupe_%a.txt
#SBATCH --array=1-5%5

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

## Locate file
repo_dir=$(git rev-parse --show-toplevel)
SAMPLE=$(awk 'BEGIN {FS="\t"} {print $1}' all_hd_samples_10_2025.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")
IMGCYT=$(awk 'BEGIN {FS="\t"} {print $2}' all_hd_samples_10_2025.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")
echo "Processing sample ${SAMPLE}"
date

## Get slide and area
SLIDE=$(echo ${SAMPLE} | cut -d "_" -f 1)
CAPTUREAREA=$(echo ${SAMPLE} | cut -d "_" -f 2)
SAM=$(paste <(echo ${SLIDE}) <(echo "-") <(echo ${CAPTUREAREA}) -d '')
echo "Slide: ${SLIDE}, capture area: ${CAPTUREAREA}"

## Find FASTQ file path
FASTQPATH=$(ls -d ${repo_dir}/raw-data/fastqs/${SAMPLE}/)

## Hank from 10x Genomics recommended setting this environment
export NUMBA_NUM_THREADS=1

spaceranger count \
    --id=${SAMPLE} \
    --transcriptome=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2024-A \
    --fastqs=${FASTQPATH} \
    --probe-set=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/Visium_Human_Transcriptome_Probe_Set_v2.1.0_GRCh38-2024-A.csv \
    --slide=${SLIDE} \
    --area=${CAPTUREAREA} \
    --cytaimage=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/raw-data/images/vis-hd/${IMGCYT}.tif \
    --image=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/raw-data/images/vis-hd/${SAMPLE}.tif \
    --create-bam=false \
    --localcores=8 \
    --localmem=64 

## Move output
echo "Moving results to new location"
date
mkdir -p ${repo_dir}/processed-data/01_spaceranger/five_samples_10_2025/
mv ${SAMPLE} ${repo_dir}/processed-data/01_spaceranger/five_samples_10_2025/

echo "**** Job ends ****"
date

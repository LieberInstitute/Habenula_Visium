#!/bin/bash
#SBATCH --mem=80G
#SBATCH -c 8
#SBATCH -p katun
#SBATCH --job-name=03-spaceranger-hd_probe_fix
#SBATCH -o logs/03-spaceranger-hd_probe_fix_%a.txt
#SBATCH -e logs/03-spaceranger-hd_probe_fix_%a.txt
#SBATCH --array=1-5%2

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## load SpaceRanger
module load spaceranger/3.1.1

## List current modules for reproducibility
module list

## Locate file
SAMPLE=$(cat 03-hd-sample-list.txt 04-hd-sample-list-250110.txt | awk 'BEGIN {FS="\t"} {print $1}' | awk "NR==${SLURM_ARRAY_TASK_ID}")
IMGCYT=$(cat 03-hd-sample-list.txt 04-hd-sample-list-250110.txt | awk 'BEGIN {FS="\t"} {print $2}' | awk "NR==${SLURM_ARRAY_TASK_ID}")

echo "Processing sample ${SAMPLE}"
date

## Get slide and area
SLIDE=$(echo ${SAMPLE} | cut -d "_" -f 1)
CAPTUREAREA=$(echo ${SAMPLE} | cut -d "_" -f 2)
SAM=$(paste <(echo ${SLIDE}) <(echo "-") <(echo ${CAPTUREAREA}) -d '')
echo "Slide: ${SLIDE}, capture area: ${CAPTUREAREA}"

## Find FASTQ file path
FASTQPATH=$(ls -d /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/raw-data/fastqs/${SAMPLE}/)

## Hank from 10x Genomics recommended setting this environment
export NUMBA_NUM_THREADS=1

spaceranger count \
    --id=${SAMPLE} \
    --transcriptome=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A \
    --fastqs=${FASTQPATH} \
    --probe-set=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A_for_lot-driven_batch_effect_correction.csv \
    --slide=${SLIDE} \
    --area=${CAPTUREAREA} \
    --cytaimage=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/raw-data/images/vis-hd/${IMGCYT}.tif \
    --image=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/raw-data/images/vis-hd/${SAMPLE}.tif \
    --loupe-alignment=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/Images/loupe-alignment/${SAM}-fiducials-image-registration.json \
    --create-bam=false \
    --localcores=8 \
    --localmem=64 


## Move output
echo "Moving results to new location"
mkdir -p /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/11_spaceranger_probe_fix
mv ${SAMPLE} /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/11_spaceranger_probe_fix/

echo "**** Job ends ****"
date

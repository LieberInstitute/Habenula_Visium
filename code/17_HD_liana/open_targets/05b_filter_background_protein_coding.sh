#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=05b_filter_background_protein_coding.R
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../../../processed-data/10_HD_bin_level/no_secondary/liana2/open_targets/logs/05b_filter_background_protein_coding.txt
#SBATCH -e ../../../../processed-data/10_HD_bin_level/no_secondary/liana2/open_targets/logs/05b_filter_background_protein_coding.txt

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load R module
module load conda_R/4.5

## List current modules for reproducibility
module list

SCRIPT_DIR="/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/code/10_HD_bin_level/liana2/open_targets"

Rscript "${SCRIPT_DIR}/05b_filter_background_protein_coding.R"

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/

#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=20G
#SBATCH --job-name=20_new_figure_LR_8combo
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../../processed-data/10_HD_bin_level/no_secondary/liana2/logs/20_new_figure_LR_8combo.txt
#SBATCH -e ../../../processed-data/10_HD_bin_level/no_secondary/liana2/logs/20_new_figure_LR_8combo.txt

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load the R module
module load conda_R/4.5

## List current modules for reproducibility
module list

## Edit with your job command
Rscript /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/code/10_HD_bin_level/liana2/20_new_figure_LR_8combo_Astros_Habenula.R
Rscript /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/code/10_HD_bin_level/liana2/20_new_figure_LR_8combo_Ependymal_MHb.R
Rscript /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/code/10_HD_bin_level/liana2/20_new_figure_LR_8combo_MHb_LHb.R
Rscript /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/code/10_HD_bin_level/liana2/20_new_figure_LR_8combo_Oligos_habenula.R

echo "**** Job ends ****"
date

#!/bin/bash
#SBATCH --partition=katun	
#SBATCH --job-name=Hb_Visium_run_all
#SBATCH -c 2
#SBATCH -t 1-00:00:00
#SBATCH --mem=30GB						                                    
#SBATCH -o logs/run_all.txt
#SBATCH -e logs/run_all.txt
# SBATCH --mail-type=ALL

# ## Explicitly pipe script output to a log
# mkdir -p logs 
log_path=logs/run_all.txt

set -e

echo "**** Job starts ****"
echo "Script RUN ALL for Habenula Visium Datasets"
echo "Samples S01 to S17"
date

echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job name (script): ${SLURM_JOB_NAME}"
echo "Hostname (computer node): ${HOSTNAME}"
echo ""
echo "Job id: ${SLURM_JOBID}"

# ## load modules
# module conda_R/4.3.x
# ## List current modules for reproducibility
# module list

MAINDIR="/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium"
CODEDIR="${MAINDIR}/code"
PROCESSEDIR="${MAINDIR}/processed-data"
PLOTDIR="${MAINDIR}/plots"


## Update code style
# cd ${CODEDIR}
# Rscript update_style.R


########  Basic workflow ########

######## Build basic spe object ########

echo "Call to process 1: build basic spe.rds object "

echo "Running 01_build_basic_spe.sh"

## Delete the logs/old-results, and re-submit EmptyDrops
cd ${CODEDIR}/02_build_spe

## Create the logs directory if it doesn't exist
[ -f logs/01_build_basic_spe.txt ] && rm logs/01_build_basic_spe.txt
# \rm logs/01_build_basic_spe.txt
mv ${PROCESSEDIR}/02_build_spe/spe.rds ${PROCESSEDIR}/02_build_spe/tmp_S1_S13

# rm ${PLOTDIR}/02_build_spe/*.png

## Build spe_raw.rds and spe.rds objects
sbatch 01_build_basic_spe.sh

echo "Call to process 1 completed!"
echo

echo "**** Job ends ****"
date



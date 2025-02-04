#!/bin/bash
#SBATCH --partition=katun	
#SBATCH --job-name=Hb_Visium_run_all
#SBATCH -c 2
#SBATCH -t 1-00:00:00
#SBATCH --mem=30GB						                                    
#SBATCH -o /dev/null
#SBATCH -e /dev/null
# SBATCH --mail-type=ALL

## Explicitly pipe script output to a log
log_path=logs/run_all.txt

{
set -e

echo "**** Job starts ****"
echo "Script RUN ALL for Habenula Visium Datasets"
echo "Samples S01 to S16"
date

echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job name (script): ${SLURM_JOB_NAME}"
echo "Hostname (computer node): ${HOSTNAME}"
echo ""
echo "Job id: ${SLURM_JOBID}"


MAINDIR="/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium"
CODEDIR="${MAINDIR}/code"
PROCESSEDIR="${MAINDIR}/processed-data"
PLOTDIR="${MAINDIR}/plots"


## Update code style
# cd ${CODEDIR}
# Rscript update_style.R


########  Basic workflow ########

######## Build basic spe object ########

echo "Running process (1) ###################################### "

## change directory
SUBDIR="02_build_spe"
cd ${CODEDIR}/${SUBDIR}

echo "Running 01_build_basic_spe.sh"

## create log dir or rm previous log files and output files
[ -f logs/01_build_basic_spe.txt ] && rm logs/01_build_basic_spe.txt

## move previous rds
[ -f ${PROCESSEDIR}/${SUBDIR}/spe.rds ] && mv ${PROCESSEDIR}/${SUBDIR}/spe.rds ${PROCESSEDIR}/${SUBDIR}/tmp
[ -f ${PROCESSEDIR}/${SUBDIR}/spe_raw.rds ] && mv ${PROCESSEDIR}/${SUBDIR}/spe_raw.rds ${PROCESSEDIR}/${SUBDIR}/tmp

echo "Previous logs and output files deleted or moved to tmp dir!"

## Run main job
id1=$(sbatch --parsable 01_build_basic_spe.sh)

echo ${id1}

echo "Process (1) completed! ###################################### "


################################################################
################################################################


echo "Running process (2) ###################################### "
echo "Running 02_scran_exploratory_QCs.sh"

## create log dir or rm previous log files and output files

[ -f logs/02_scran_exploratory_QCs.txt ] && rm logs/02_scran_exploratory_QCs.txt
rm ${PLOTDIR}/${SUBDIR}/*.pdf 

## Run dependency job
id2=$(sbatch --parsable --dependency=afterok:$id1 02_scran_exploratory_QCs.sh)

echo $id2

echo "Process (2) completed! ###################################### "


echo "**** Job ends ****"
date

} > $log_path 2>&1

## Cynthia SC - Feb, 2025

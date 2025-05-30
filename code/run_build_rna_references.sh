#!/bin/bash
#SBATCH --partition=katun
#SBATCH --job-name=run_build_rna_references
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH --mem=15GB
#SBATCH -o /dev/null
#SBATCH -e /dev/null
# SBATCH --mail-type=ALL

## Explicitly pipe script output to a log
log_path=logs/run_build_rna_references.txt

{
set -e

echo "**** Job starts ****"
echo "Build single-cell RNAseq and multiome references"
date

echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job name (script): ${SLURM_JOB_NAME}"
echo "Hostname (computer node): ${HOSTNAME}"
echo ""
echo "Job id: ${SLURM_JOBID}"



echo "Main Directories ########################################################## "

MAINDIR="/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium"
CODEDIR="${MAINDIR}/code"
PROCESSEDIR="${MAINDIR}/processed-data"
PLOTDIR="${MAINDIR}/plots"

echo "Main dir: ${MAINDIR}"
echo "Processed dir: ${PROCESSEDIR}"
echo "Plot dir: ${PLOTDIR}"


## Build enrichment stats objects fro
SUBDIR="05_snRNA-seq_model_stats"
cd ${CODEDIR}/${SUBDIR}
echo "Current code dir: ${CODEDIR}/${SUBDIR}/"

echo "Build snRNASeq reference with final annotations from: "
echo "https://github.com/LieberInstitute/Habenula_Pilot ####################### "
rm -f logs/01_pseudobulk_reference.txt
rm -f ${PROCESSEDIR}/${SUBDIR}/pseudobulk_final_*.rds
rm -f ${PROCESSEDIR}/${SUBDIR}/enrichment_final_*.rds
Rscript 01_pseudobulk_reference.R

echo "Build multiome-RNA reference with final annotations from: "
echo "https://github.com/LieberInstitute/Hb_multiome ######################### "
rm -f logs/02_multiome_rna_reference.txt
#rm -f ${PROCESSEDIR}/${SUBDIR}/
rm -f ${PROCESSEDIR}/${SUBDIR}/enrichment_snRNA-multiome*.rds
sbatch 02_multiome_rna_reference.sh

echo "Build multiome-RNA reference with final annotations from: "
echo "https://github.com/LieberInstitute/Hb_multiome ######################### "
echo "Alternative version for Habenula clusters merged"
rm -f logs/03_pseudobulk_reference_habenula_merged.txt
rm -f ${PROCESSEDIR}/${SUBDIR}/pseudobulk_Hb_merged_final*.rds
rm -f ${PROCESSEDIR}/${SUBDIR}/enrichment_Hb_merged_final*.rds
sbatch 03_pseudobulk_reference_habenula_merged.sh


echo "**** Job ends ****"
date

} > $log_path 2>&1

## Cynthia SC - May, 2025

#!/bin/bash
#SBATCH --partition=katun
#SBATCH --job-name=run_build_rna_references
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH --mem=2GB
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


## Build enrichment stats objects
SUBDIR="05_snRNA-seq_model_stats"
cd ${CODEDIR}/${SUBDIR}
echo "Current code dir: ${CODEDIR}/${SUBDIR}/"

echo "Build snRNASeq reference with final annotations from: "
echo "https://github.com/LieberInstitute/Habenula_Pilot ####################### "
mv logs/01_pseudobulk_reference.txt logs/old/ 2>/dev/null || true
rm -f ${PROCESSEDIR}/${SUBDIR}/pseudobulk_final_*.rds
rm -f ${PROCESSEDIR}/${SUBDIR}/enrichment_final_*.rds

id=$(sbatch --parsable 01_pseudobulk_reference.sh)
#ls -1t "${PROCESSEDIR}/${SUBDIR}/"
sbatch --dependency=afterok:$id --wrap="ls -1t \"${PROCESSEDIR}/${SUBDIR}/\""

echo "Build multiome-RNA reference with final annotations from: "
echo "Multiome reference is run on Multiome project"


echo "**** Job ends ****"
date

} > $log_path 2>&1

## Cynthia SC - Aug, 2025

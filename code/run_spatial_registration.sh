#!/bin/bash
#SBATCH --partition=katun
#SBATCH --job-name=run_spatial_registration
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH --mem=15GB
#SBATCH -o /dev/null
#SBATCH -e /dev/null
# SBATCH --mail-type=ALL

## Explicitly pipe script output to a log
log_path=logs/run_spatial_registration.txt

{
set -e

echo "**** Job starts ****"
echo "Run Spatial-Registration"
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



echo "#########   Spatial Registrattion              ############################# "
echo "#########   - snRNAseq vs Visium               ############################# "
echo "#########   - snRNAseq vs Multiome RNA.        ############################# "


## change directory
SUBDIR="05_brain_area_differential_expression"
cd ${CODEDIR}/${SUBDIR}
echo "Current code dir: ${CODEDIR}/${SUBDIR}/"

## we do not need to run the clustering again if only reference change, so you can safely skip this process
# 
# ##### Compute pseudobulk and correlation model for Visium BayesSpace k=(2-28)
# ## The model statistics are required to compute spatial-registration, if done, jump this chunck ------------
# 
# rm -f logs/01_create_pseudobulk_data_*.txt 
# rm -f ${PROCESSEDIR}/${SUBDIR}/stats_summary_csv/*_basic_stats.csv
# rm -f ${PROCESSEDIR}/${SUBDIR}/stats_summary_csv/*_SpatialD_info.csv
# rm -f ${PROCESSEDIR}/${SUBDIR}/sce_pseudo_PCA_brain_area*.rds
# 
# id1=$(sbatch --parsable 01_create_pseudobulk_data.sh)
# echo "Running array job: ${id1}"
# 
# # dependency array job 
# rm -f logs/06_model_BayesSpace_*.txt
# rm -f ${PROCESSEDIR}/${SUBDIR}/modeling_results_BS/modeling_results_BayesSpace*.Rdata
# id2=$(sbatch --parsable --dependency=afterok:$id1 06_model_BayesSpace.sh)
# echo "Running dependency array job: ${id2}"
# 
# #------------------------------------------------------------------------------------------------------------
# 
# ##### EDA: build HISTOGRAMS and STACKED BAR PLOTS with Hb and no-Habenula count/proportions by BayesSpace domain
# # Here we compare Hb-Taxomony manual annotations (RNAScope) vs SpD in clustering
# 
# # independent job
# rm -f logs/08_manual_ann_vs_bayes_space_*.txt
# rm -f ${PLOTDIR}/${SUBDIR}/08_manual_ann_vs_bayes_space/*.pdf
# sbatch 08_manual_ann_vs_bayes_space.sh
# 
# #------------------------------------------------------------------------------------------------------------

# submit dummy Slurm jobs that do nothing but exit successfully to be able to skip the lines above
id1=$(sbatch --parsable --wrap="echo 'Skipping 06_model_BayesSpace.sh'; sleep 1")
id2=$(sbatch --parsable --dependency=afterok:$id1 \
    --wrap="echo 'Skipping 06_model_BayesSpace.sh'; sleep 1")
echo "Dummy job submitted with ID: ${id2}"


#####  Compute Spatial registration for VISIUM Bayes-Space vs both FINE and BROAD (snRNAseq)
## x-axis = snRNAseq cell-types
## y-axis = spatial Habenula Visium domains

## change directory
SUBDIR="06_spatial_registration_vs_snRNA-seq"
cd ${CODEDIR}/${SUBDIR}
echo "Current code dir: ${CODEDIR}/${SUBDIR}/"

# dependency array job 
echo " Spatial Registrattion Visum  vs scRNAseq human pilot"

mv logs/01_compute_cor.* logs/old/ 2>/dev/null || true
if [ -f "${PROCESSEDIR}/${SUBDIR}/cor_BayesSpace_vs_snRNA-seq_top100.Rdata" ]; then
    mkdir -p old
    mv "${PROCESSEDIR}/${SUBDIR}/cor_BayesSpace_vs_snRNA-seq_top100.Rdata" old/
fi
if [ -f "${PLOTDIR}/${SUBDIR}/*_broadRes.pdf" ]; then
    mkdir -p old
    mv "${PLOTDIR}/${SUBDIR}/*_broadRes.pdf" old/
fi
if [ -f "${PLOTDIR}/${SUBDIR}/*_fineRes.pdf" ]; then
    mkdir -p old
    mv "${PLOTDIR}/${SUBDIR}/*_fineRes.pdf" old/
fi

id3=$(sbatch --parsable --dependency=afterok:$id2 01_compute_cor.sh)
echo "Running dependency array job: ${id3}"


echo " Spatial Registrattion Visium vs scRNAseq human pilot with Hb clusters merged"
## First remove old data and plots
rm -f logs/02_compute_corr_hb_merged.*.out
rm -f logs/02_compute_corr_hb_merged.*.err
rm -f ${PROCESSEDIR}/${SUBDIR}/cor_BayesSpace_vs_snRNA-seq_top100_Hb_merged.Rdata
rm -f ${PLOTDIR}/${SUBDIR}/*_broadRes_Hb_merged.pdf

sbatch 02_compute_corr_hb_merged.sh


#------------------------------------------------------------------------------------------------------------

echo " Spatial Registrattion scRNAseq human pilot vs multiome-RNA human hb"

SUBDIR="07_spatial_registration_vs_multiome_snRNA-seq"
cd ${CODEDIR}/${SUBDIR}
echo "Current code dir: ${CODEDIR}/${SUBDIR}/"

## Compute correlations for visium vs snRNAseq fine and broad resolution (CSC) 
# Plot in both verical and horizontal formatR
rm -f logs/01_compute_cor_snRnaseq_multiomeRnaseq_*.txt
rm -f ${PROCESSEDIR}/${SUBDIR}/cor_multiome_vs_snRNA-seq_top100_*.Rdata
rm -f ${PLOTDIR}/${SUBDIR}/cor_top100_registration_snMultiome_snRNAseq_v2_*.pdf

sbatch 01_compute_cor_snRnaseq_multiomeRnaseq.sh

## Compute correlations for visium vs snRNAseq fine and broad resolution (CSC)
rm -f logs/02_compute_cor_visium_multiomeRnaseq.txt
rm -f ${PROCESSEDIR}/${SUBDIR}/bayesSpace_cor_top100_*.Rdata
rm -f ${PLOTDIR}/${SUBDIR}/cor_top100_spatial_registration_snMultiome_v2.pdf

sbatch 02_compute_cor_visium_multiomeRnaseq.sh

echo "**** Job ends ****"
date

} > $log_path 2>&1

## Cynthia SC - May, 2025

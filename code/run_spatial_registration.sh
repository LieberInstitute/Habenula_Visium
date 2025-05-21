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
echo "Build single-cell RNA-seq and multiome-RNA references for Spatial-Registration"
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


##### Compute pseudobulk and correlation model for Visium BayesSpace k=(2-28)
## The model statistics are required to compute spatial-registration, if done, jump this chunck ------------

rm -f logs/01_create_pseudobulk_data_*.txt
rm -f ${PROCESSEDIR}/${SUBDIR}/stats_summary_csv/*_basic_stats.csv
rm -f ${PROCESSEDIR}/${SUBDIR}/stats_summary_csv/*_SpatialD_info.csv
rm -f ${PROCESSEDIR}/${SUBDIR}/sce_pseudo_PCA_brain_area*.rds
id1=$(sbatch --parsable 01_create_pseudobulk_data.sh)
echo "Running array job: ${id1}"

# dependency job 
rm -f logs/06_model_BayesSpace_*.txt
rm -f ${PROCESSEDIR}/${SUBDIR}/modeling_results_BS/modeling_results_BayesSpace*.Rdata
#sbatch 06_model_BayesSpace.sh
id2=$(sbatch --parsable --dependency=afterok:$id1 06_model_BayesSpace.sh)
echo "Running dependency array job: ${id2}"

#------------------------------------------------------------------------------------------------------------

##### EDA: build HISTOGRANS and STACKED BAR PLOTS with Hb and no-Habenula count/proportions by BayesSpace domain
# Here we compare Hb-Taxomony manual annotations (RNAScope) vs SpD in clustering

rm -f logs/08_manual_ann_vs_bayes_space_*.txt
rm -f ${PLOTDIR}/${SUBDIR}/08_manual_ann_vs_bayes_space/*.pdf
sbatch 08_manual_ann_vs_bayes_space.sh

#------------------------------------------------------------------------------------------------------------

#####  Compute Spatial registration for VISIUM Bayes-Space vs both FINE and BROAD (snRNAseq)
## x-axis = snRNAseq cell-types
## y-axis = spatial Habenula Visium domains

## change directory
SUBDIR="06_spatial_registration_vs_snRNA-seq"
cd ${CODEDIR}/${SUBDIR}
echo "Current code dir: ${CODEDIR}/${SUBDIR}/"


echo " Spatial Registrattion Visum  vs scRNAseq human pilot"
rm -f logs/01_compute_cor.*.out
rm -f logs/01_compute_cor.*.err
rm -f ${PROCESSEDIR}/${SUBDIR}/cor_BayesSpace_vs_snRNA-seq_top100.Rdata
rm -f ${PLOTDIR}/${SUBDIR}/*_broadRes.pdf
rm -f ${PLOTDIR}/${SUBDIR}/*_fineRes.pdf
sbatch 01_compute_cor.sh


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

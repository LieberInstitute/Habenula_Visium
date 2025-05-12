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


echo "Main Directories ########################################################## "

MAINDIR="/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium"
CODEDIR="${MAINDIR}/code"
PROCESSEDIR="${MAINDIR}/processed-data"
PLOTDIR="${MAINDIR}/plots"

echo "Main dir: ${MAINDIR}"
echo "Processed dir: ${PROCESSEDIR}"
echo "Plot dir: ${PLOTDIR}"



echo "*** NOTES: ***"
echo "To re-build the single-cell references go to:  `run_build_rna_references.sh`"
echo "To re-build spatial-registration go to: `run_spatial_registration.sh`"



echo "Build spe basic ########################################################## "

SUBDIR="02_build_spe"
cd ${CODEDIR}/${SUBDIR}
echo "Current code dir: ${CODEDIR}/${SUBDIR}/"

## create log dir or rm previous log files and output files
[ -f logs/01_build_basic_spe.txt ] && rm logs/01_build_basic_spe.txt
## move previous rds
[ -f ${PROCESSEDIR}/${SUBDIR}/spe.rds ] && mv ${PROCESSEDIR}/${SUBDIR}/spe.rds ${PROCESSEDIR}/${SUBDIR}/tmp
[ -f ${PROCESSEDIR}/${SUBDIR}/spe_raw.rds ] && mv ${PROCESSEDIR}/${SUBDIR}/spe_raw.rds ${PROCESSEDIR}/${SUBDIR}/tmp
echo "Previous logs and output files deleted or moved to tmp dir!"

id1=$(sbatch --parsable 01_build_basic_spe.sh)
echo ${id1}


echo "Scran exploratory QCs #################################################### "

[ -f logs/02_scran_exploratory_QCs.txt ] && rm logs/02_scran_exploratory_QCs.txt
rm -f ${PLOTDIR}/${SUBDIR}/*.pdf
## mv previous rds
[ -f ${PROCESSEDIR}/${SUBDIR}/spe_qc_low_lib_edge.rds ] && rm ${PROCESSEDIR}/${SUBDIR}/spe_qc_low_lib_edge.rds

## Run dependency job
id2=$(sbatch --parsable --dependency=afterok:$id1 02_scran_exploratory_QCs.sh)
echo $id2


echo "SpotSweeper QCs ######################################################### "

[ -f logs/03_SpotSweeper.txt ] && rm logs/03_SpotSweeper.txt
rm -f ${PLOTDIR}/${SUBDIR}/03_SpotSweeper/*.pdf
rm -f ${PROCESSEDIR}/${SUBDIR}/03_SpotSweeper/*.csv
## mv previous rds
[ -f ${PROCESSEDIR}/${SUBDIR}/spe_scran_spotsweeper.rds ] && rm ${PROCESSEDIR}/${SUBDIR}/spe_scran_spotsweeper.rds

## Run dependency job after scran outliers identification
id3=$(sbatch --parsable --dependency=afterok:$id2 03_SpotSweeper.sh)
echo $id3


echo "#########   Batch correction and clustering    ########################## "

######## Prepare dataset, normalize, harmonize and run BayesSpace ##############

SUBDIR="04_harmony_BayesSpace"
cd ${CODEDIR}/${SUBDIR}
echo "Current code dir: ${CODEDIR}/${SUBDIR}/"

01-filter_normalize.R
02-compute_GLM-PCA.sh

## Run Harmony

## First remove old plots
rm -f logs/03-preprocess_and_harmony.txt
rm -f ${PLOTDIR}/${SUBDIR}/reduction_dimension_PCA.pdf
rm -f ${PLOTDIR}/${SUBDIR}/reduction_dimension_GLMPCA.pdf
rm -f ${PLOTDIR}/${SUBDIR}/reduction_dimension_GLMPCA_other_features.pdf
rm -f ${PLOTDIR}/${SUBDIR}/tSNE_perplexity*.pdf
rm -f ${PLOTDIR}/${SUBDIR}/UMAP_*.pdf
rm -f ${PLOTDIR}/${SUBDIR}/reduction_dimension_Harmony_vs_GLMPCA.pdf
rm -f ${PLOTDIR}/${SUBDIR}/reduction_dimension_Harmony_vs_GLMPCA_sumUMI_sumGene.pdf
rm -f ${PLOTDIR}/${SUBDIR}/graph_based_harmony.pdf
## and, remove old data objects
rm -rf ${PROCESSEDIR}/${SUBDIR}/clusters_graphbased/
rm -rf ${PROCESSEDIR}/${SUBDIR}/clusters_graphbased_cut_at/
rm -rf ${PROCESSEDIR}/${SUBDIR}/g*harmony.Rdata
rm -rf ${PROCESSEDIR}/${SUBDIR}/spe_harmony.rds

## Run job
sbatch 03-preprocess_and_harmony.sh # need to add the sh script


echo "#########   BayesSpace ################################################## "

## create log dir or rm previous log files and output files
rm -f logs/04-BayesSpace_k_search*.txt
rm -f ${PLOTDIR}/${SUBDIR}/BayesSpace/BayesSpace_harmony_k*_raw.pdf
rm -rf ${PROCESSEDIR}/${SUBDIR}/clusters_BayesSpace/BayesSpace_harmony_k*

sbatch 04-BayesSpace_k_search.sh

echo "#########   cell-type differential expression  ########################### "

## change directory
SUBDIR="05_layer_differential_expression"
cd ${CODEDIR}/${SUBDIR}
echo "Current code dir: ${CODEDIR}/${SUBDIR}/"

## First remove old plots
rm -f logs/01_create_pseudobulk_data_*.txt
# rm -f ${PLOTDIR}/${SUBDIR}/sce_pseudo_gene_explanatory_vars_k*.pdf
rm -f ${PROCESSEDIR}/${SUBDIR}/sce_pseudo_BayesSpace_k*.rds

sbatch 01_create_pseudobulk_data.sh


echo "#########   Running exploring variance ################################## "

## First remove old plots
rm -f logs/02_explore_expr_variability_*.txt
rm -f ${PLOTDIR}/${SUBDIR}/sce_pseudo_PCs_k*.pdf
rm -f ${PLOTDIR}/${SUBDIR}/sce_pseudo_gene_explanatory_vars_k*.pdf

sbatch 02_explore_expr_variability.sh


echo  "#########   Running model BayesSpace #################################### "

## First remove old plots
rm -f logs/03_model_BayesSpace*.err
rm -f logs/03_model_BayesSpace*.out

sbatch 03_model_BayesSpace.sh


echo "#########   Explore BRAIN-AREA differential expression  ################# "

## change directory
SUBDIR="05_brain_area_differential_expression"
cd ${CODEDIR}/${SUBDIR}
echo "Current code dir: ${CODEDIR}/${SUBDIR}/"

## First remove old plots
rm -f logs/01_create_pseudobulk_data_BS*.txt
rm -f ${PROCESSEDIR}/${SUBDIR}/sce_pseudo_PCA_brain_area_k*.rds
rm -f ${PROCESSEDIR}/${SUBDIR}/stats_summary_csv/*.csv

id1_pseudoDE=$(sbatch --parsable 01_create_pseudobulk_data.sh)

## First remove old plots
# Note: Previous results using the snRNA-seq reference with Hb split into MHb and LHb are stored in old_hb_no-merged/.
rm -f logs/02_explore_expr_variability_BS*.txt
rm -f ${PLOTDIR}/${SUBDIR}/sce_pseudo_PC*.pdf
rm -f ${PLOTDIR}/${SUBDIR}/sce_pseudo_gene_explanatory_vars_k*.pdf

#sbatch 02_explore_expr_variability.sh
id2_pseudoDE=$(sbatch --parsable --dependency=afterok:$id1_pseudoDE 02_explore_expr_variability.sh)

## First remove old plots
rm -f logs/03_covariate_analysis.txt
rm -f ${PLOTDIR}/${SUBDIR}/03_covariate_analysis/*.pdf
rm -f ${PLOTDIR}/${SUBDIR}/03_covariate_analysis/*.png

#sbatch 02_explore_expr_variability.sh
# sbatch 03_covariate_analysis.sh
sbatch --dependency=afterok:$id2_pseudoDE 03_covariate_analysis.sh



echo "#########   Make SpatialLIBD app.              ############################# "
echo "#########   For pseudobulk data.               ############################# "

## (1) Be sure to create a smaller spe object for the shiny app
## Ex. spe_harmony_shiny.rds or spe_pseudobulk_shiny.rds
code/03_spatialLIBD_app_pseudobulk/03_spatialLIBD_app_pseudobulk/01_make_spe_shiny_app.R
## (2) create soft-links, dirs, readme and documentation
code/03_spatialLIBD_app_pseudobulk/03_spatialLIBD_app_pseudobulk/select_BayesSpaceK_brainModel.R
## (3) check/test app
code/03_spatialLIBD_app_pseudobulk/03_spatialLIBD_app_pseudobulk/app.R
## (3) check/test app
code/03_spatialLIBD_app_pseudobulk/03_spatialLIBD_app_pseudobulk/deploy.R

echo "**** Job ends ****"
date

} > $log_path 2>&1

## Cynthia SC - Feb, 2025

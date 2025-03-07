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

echo "Main dir: ${MAINDIR}"
echo "Processed dir: ${PROCESSEDIR}"
echo "Plot dir: ${PLOTDIR}"


## Update code style
# cd ${CODEDIR}
# Rscript update_style.R


########  Basic workflow ########

######## Build basic spe object ########

echo "Running process (1) ###################################### "

## change directory
SUBDIR="02_build_spe"

cd ${CODEDIR}/${SUBDIR}

echo "Current code dir: ${CODEDIR}/${SUBDIR}/"

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

echo "Build spe completed! ###################################### "

echo "Running process (2) ###################################### "
echo "Running 02_scran_exploratory_QCs.sh"

## create log dir or rm previous log files and output files

[ -f logs/02_scran_exploratory_QCs.txt ] && rm logs/02_scran_exploratory_QCs.txt
rm -f ${PLOTDIR}/${SUBDIR}/*.pdf

## mv previous rds
[ -f ${PROCESSEDIR}/${SUBDIR}/spe_qc_low_lib_edge.rds ] && rm ${PROCESSEDIR}/${SUBDIR}/spe_qc_low_lib_edge.rds

## Run dependency job
id2=$(sbatch --parsable --dependency=afterok:$id1 02_scran_exploratory_QCs.sh)

echo $id2

echo "EDA with scran completed! ###################################### "


echo "Running process (3) ###################################### "
echo "Running 03_SpotSweeper.sh"

## create log dir or rm previous log files and output files

[ -f logs/03_SpotSweeper.txt ] && rm logs/03_SpotSweeper.txt
rm ${PLOTDIR}/${SUBDIR}/03_SpotSweeper/*.pdf
#[ -f ${PROCESSEDIR}/${SUBDIR}/03_SpotSweeper/*.csv ] &&
rm ${PROCESSEDIR}/${SUBDIR}/03_SpotSweeper/*.csv

## mv previous rds
[ -f ${PROCESSEDIR}/${SUBDIR}/spe_scran_spotsweeper.rds ] && rm ${PROCESSEDIR}/${SUBDIR}/spe_scran_spotsweeper.rds

## Run dependency job after scran outliers identification
id3=$(sbatch --parsable --dependency=afterok:$id2 03_SpotSweeper.sh)

echo $id3

echo "EDA with SpotSweeper completed! ###################################### "


echo "**** Job ends ****"
date

} > $log_path 2>&1


echo "############################################################################ "
echo "#########   Batch correction and clustering    ############################# "
echo "#########                                      ############################# "
echo "############################################################################ "

echo "Running Harmoy ###################################### "
echo "sbatch 03-preprocess_and_harmony.sh"

######## Prepare dataset, normalize, harmonize and run BayesSpace ########

## change directory
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

echo "Previous logs and output files deleted!"

## Run job
sbatch 03-preprocess_and_harmony.sh # need to add the sh script

echo "Processing Harmony ... wait"
squeue -u csoto

echo "Harmony completed! ###################################### "


echo "Running BayesSpace ###################################### "
echo "04-BayesSpace_k_search.sh"

## create log dir or rm previous log files and output files
rm -f logs/04-BayesSpace_k_search*.txt
rm -f ${PLOTDIR}/${SUBDIR}/BayesSpace/BayesSpace_harmony_k*_raw.pdf
rm -rf ${PROCESSEDIR}/${SUBDIR}/clusters_BayesSpace/BayesSpace_harmony_k*

echo "Previous logs and output files deleted!"

## Run job
sbatch 04-BayesSpace_k_search.sh

echo "Processing BayesSpace ... wait"
squeue -u csoto

echo "BayesSpace ###################################### "


echo "############################################################################ "
echo "#########   Layer differential expression.     ############################# "
echo "#########                                      ############################# "
echo "############################################################################ "

echo "Running pseudobulk ###################################### "
echo "01_create_pseudobulk_data.sh"

## change directory
SUBDIR="05_layer_differential_expression"

cd ${CODEDIR}/${SUBDIR}

echo "Current code dir: ${CODEDIR}/${SUBDIR}/"

## First remove old plots
rm -f logs/01_create_pseudobulk_data_*.txt
rm -f ${PLOTDIR}/${SUBDIR}/sce_pseudo_gene_explanatory_vars_k*.pdf
rm -f ${PROCESSEDIR}/${SUBDIR}/sce_pseudo_BayesSpace_k*.rds

sbatch 01_create_pseudobulk_data.sh
squeue -u csoto


echo "Running exploring variance ###################################### "
echo "02_explore_expr_variability.sh"

## First remove old plots
rm -f logs/02_explore_expr_variability_*.txt
rm -f ${PLOTDIR}/${SUBDIR}/sce_pseudo_PCs_k*.pdf
rm -f ${PLOTDIR}/${SUBDIR}/sce_pseudo_gene_explanatory_vars_k*.pdf

sbatch 02_explore_expr_variability.sh
squeue -u csoto


echo "Running model BayesSpace ###################################### "
echo "03_model_BayesSpace.sh"

## First remove old plots
rm -f logs/03_model_BayesSpace*.err
rm -f logs/03_model_BayesSpace*.out

sbatch 03_model_BayesSpace.sh



echo "############################################################################ "
echo "#########   Compute enrichment with registration_wrapper  ################## "
echo "#########   From RNA multiome modality (WNN Leiden res=2, knn=30) ########## "
echo "############################################################################ "

## change directory
SUBDIR="05_snRNA-seq_model_stats"

cd ${CODEDIR}/${SUBDIR}

echo "Current code dir: ${CODEDIR}/${SUBDIR}/"

## First remove old plots
rm -f logs/02_multiome_rna_reference.txt
rm -f ${PROCESSEDIR}/${SUBDIR}/enrichment_snRNA-multiome.rds

sbatch 02_multiome_rna_reference.sh



echo "############################################################################ "
echo "#########   Spatial Registrattion              ############################# "
echo "#########                                      ############################# "
echo "############################################################################ "


## change directory
SUBDIR="06_spatial_registration_vs_snRNA-seq"

cd ${CODEDIR}/${SUBDIR}

echo "Current code dir: ${CODEDIR}/${SUBDIR}/"

## First remove old plots
rm -f logs/01_compute_cor.*.out
rm -f logs/01_compute_cor.*.err
rm -f ${PROCESSEDIR}/${SUBDIR}/cor_BayesSpace_vs_snRNA-seq.Rdata

sbatch 01_compute_cor.sh



echo "Compute layer correlation annotation ###################################### "

echo "Current code dir: ${CODEDIR}/${SUBDIR}/"

rm -f ${PROCESSEDIR}/${SUBDIR}/bayesSpacce_layer_cor_top100.Rdata
rm -f ${PLOTDIR}/${SUBDIR}/cor_top100_spatial_registration.pdf

Rscript 01_layer_correlation_annotation.R



echo "Plot correlation  ######################################################### "

echo "Current code dir: ${CODEDIR}/${SUBDIR}/"

rm -f ${PLOTDIR}/${SUBDIR}/snRNA-seq_registration_fineRes_*.pdf
rm -f ${PLOTDIR}/${SUBDIR}/snRNA-seq_registration_broadRes_*.pdf

Rscript 02_plot_cor_basic.R







## Cynthia SC - Feb, 2025

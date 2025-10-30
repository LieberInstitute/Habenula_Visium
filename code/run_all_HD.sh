#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=2G
#SBATCH --job-name=run_all_HD
#SBATCH -c 1
#SBATCH -t 10:00
#SBATCH -o ../run_all_HD.txt
#SBATCH -e ../run_all_HD.txt

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

repo_dir=$(git rev-parse --show-toplevel)

################################################################################
#   Create SpatialExperiments
################################################################################

#   Bin-level
cd $repo_dir/code/10_HD_bin_level
job_id_1_1=$(sbatch --parsable 01_build_spe.sh)
job_id_1_2=$(sbatch --dependency=afterok:${job_id_1_1} --parsable 02_QC.sh)

#   Cell-level
cd $repo_dir/code/09_HD_cell_level
job_id_1_3=$(sbatch --parsable 01_bin2cell.sh)
job_id_1_4=$(sbatch --dependency=afterok:${job_id_1_3} --parsable 02_build_spe_raw.sh)
job_id_1_5=$(sbatch --dependency=afterok:${job_id_1_4} --parsable 03_build_spe_QC.sh)
job_id_1_6=$(sbatch --dependency=afterok:${job_id_1_5} --parsable 03_new_sample_report.R)
job_id_1_7=$(sbatch --dependency=afterok:${job_id_1_6} --parsable 04_HVG.sh)

################################################################################
#   FICTURE with library-size-normalized inputs
################################################################################

cd $repo_dir/code/10_HD_bin_level/ficture_harmony

#-------------------------------------------------------------------------------
#   Main steps to run FICTURE
#-------------------------------------------------------------------------------

job_id_2_1=$(sbatch --parsable 01_build_spe.sh)
job_id_2_2=$(sbatch --dependency=afterok:${job_id_2_1} --parsable 02_normalized_input.sh)
job_id_2_3=$(sbatch --dependency=afterok:${job_id_2_2} --parsable 03_ficture_run.sh)
job_id_2_4=$(sbatch --dependency=afterok:${job_id_2_3} --parsable 03_ficture_plot.sh)
job_id_2_5=$(sbatch --dependency=afterok:${job_id_2_4} --parsable 04_spatula_join.sh)
job_id_2_6=$(sbatch --dependency=afterok:${job_id_2_5} --parsable 09_bin_level_merge_norm.sh)

#-------------------------------------------------------------------------------
#   Spatial registration
#-------------------------------------------------------------------------------

job_id_2_7=$(sbatch --dependency=afterok:${job_id_2_6} --parsable 10_registration_wrapper_norm.sh)
job_id_2_8=$(sbatch --dependency=afterok:${job_id_2_7} --parsable 11_cor_heatmap_norm.sh)

################################################################################
#   FICTURE with cleaningY batch-corrected inputs
################################################################################

#-------------------------------------------------------------------------------
#   Main steps to run FICTURE
#-------------------------------------------------------------------------------

job_id_3_1=$(sbatch --dependency=afterok:${job_id_2_1} --parsable 08_batcheffect_lm.sh)
job_id_3_2=$(sbatch --dependency=afterok:${job_id_3_1} --parsable 08_batcheffect_com.sh)
job_id_3_3=$(sbatch --dependency=afterok:${job_id_3_2} --parsable 08_reprepareinput_ficture.sh)
job_id_3_4=$(sbatch --dependency=afterok:${job_id_3_3} --parsable 08_rerun_ficture.sh)
job_id_3_5=$(sbatch --dependency=afterok:${job_id_3_4} --parsable 03_ficture_plot_cleany.sh)
job_id_3_6=$(sbatch --dependency=afterok:${job_id_3_5} --parsable 09_get_ficture_clusters.sh)
job_id_3_7=$(sbatch --dependency=afterok:${job_id_3_6} --parsable 09_bin_level_merge_cleany.sh)

#-------------------------------------------------------------------------------
#   Spatial registration
#-------------------------------------------------------------------------------

job_id_3_8=$(sbatch --dependency=afterok:${job_id_3_7} --parsable 10_registration_wrapper_cleany.sh)
job_id_3_9=$(sbatch --dependency=afterok:${job_id_3_8} --parsable 11_cor_heatmap_cleany.sh)

################################################################################
#   Banksy
################################################################################

#   Finding SVGs
cd $repo_dir/code/10_HD_bin_level
job_id_4_1=$(sbatch --dependency=afterok:${job_id_1_6} --parsable 03_rasterize.sh)
job_id_4_2=$(sbatch --dependency=afterok:${job_id_4_1} --parsable 04_nnSVG.sh)
job_id_4_3=$(sbatch --dependency=afterok:${job_id_4_2} --parsable 05_gather_variable_genes.sh)

#   Running Banksy
cd $repo_dir/code/09_HD_cell_level
job_id_4_4=$(sbatch --dependency=afterok:${job_id_4_3} --parsable 05_banksy_embedding.sh)
job_id_4_5=$(sbatch --dependency=afterok:${job_id_4_4} --parsable 06_banksy_clustering.sh)

#-------------------------------------------------------------------------------
#   Spatial registration
#-------------------------------------------------------------------------------

cd $repo_dir/code/09_HD_cell_level/registration_banksy
job_id_4_6=$(sbatch --dependency=afterok:${job_id_4_5} --parsable 01_registration_wrapper.sh)
job_id_4_7=$(sbatch --dependency=afterok:${job_id_4_6} --parsable 02_cor_heatmap.sh)
job_id_4_8=$(sbatch --dependency=afterok:${job_id_4_7} --parsable 06_deciding_k.sh)


################################################################################
#   Other cell-level analyses
################################################################################

cd $repo_dir/code/09_HD_cell_level
job_id_5_1=$(sbatch --dependency=afterok:${job_id_4_8} --parsable 07_xenium_genes.sh)

cd $repo_dir/code/09_HD_cell_level/quick_shiny
#   Now locally run 01_app.R, exporting CSV of habenula and thalamus annotations

cd $repo_dir/code/09_HD_cell_level
job_id_5_2=$(sbatch --dependency=afterok:${job_id_4_8} --parsable 14_crawdad_region_prep.sh)
job_id_5_3=$(sbatch --dependency=afterok:${job_id_5_2} --parsable 15_crawdad_region_run.sh)
job_id_5_4=$(sbatch --dependency=afterok:${job_id_5_3} --parsable 16_crawdad_region_plot.sh)

################################################################################
#   Other bin-level analyses
################################################################################

cd $repo_dir/code/10_HD_bin_level/ficture_harmony
job_id_6_1=$(sbatch --dependency=afterok:${job_id_4_8} --parsable 12_sample_specificity.sh)
job_id_6_2=$(sbatch --dependency=afterok:${job_id_6_1} --parsable 13_jaccard.sh)
job_id_6_3=$(sbatch --dependency=afterok:${job_id_6_2} --parsable 14_batch_sum.sh)
job_id_6_4=$(sbatch --dependency=afterok:${job_id_6_3} --parsable 15_libsize_sum.sh)
job_id_6_5=$(sbatch --dependency=afterok:${job_id_6_4} --parsable 16_libsize_batch.sh)
job_id_6_6=$(sbatch --dependency=afterok:${job_id_6_5} --parsable 17_habenula_spot_plots.sh)
job_id_6_7=$(sbatch --dependency=afterok:${job_id_6_6} --parsable 18_sum_plot.sh)
job_id_6_8=$(sbatch --dependency=afterok:${job_id_6_7} --parsable 19_4_metrics_barplot.sh)
job_id_6_9=$(sbatch --dependency=afterok:${job_id_6_8} --parsable 19_4_metrics_lineplot.sh)
job_id_6_10=$(sbatch --dependency=afterok:${job_id_6_9} --parsable 20_plot_paired_data.sh)

cd $repo_dir/code/10_HD_bin_level
job_id_6_11=$(sbatch --dependency=afterok:${job_id_1_6} --parsable 06_plot_markers.sh)

################################################################################
#   LIANA
################################################################################

cd $repo_dir/code/10_HD_bin_level/exploratory
job_id_7_1=$(sbatch --dependency=afterok:${job_id_4_8} --parsable 12_preprocess_anndata.sh)
job_id_7_2=$(sbatch --dependency=afterok:${job_id_7_1} --parsable 14_0_liana_prepare_inputs.sh)
job_id_7_3=$(sbatch --dependency=afterok:${job_id_7_2} --parsable 14_2_LIANA+_preprocess.sh)
job_id_7_4=$(sbatch --dependency=afterok:${job_id_7_3} --parsable 14_3_LIANA+_preprocss_extracellular.sh)
job_id_7_5=$(sbatch --dependency=afterok:${job_id_7_4} --parsable 15_1_LIANA+_cell.sh)
job_id_7_6=$(sbatch --dependency=afterok:${job_id_7_4} --parsable 15_1_LIANA+_extracell.sh)
job_id_7_7=$(sbatch --dependency=afterok:${job_id_7_5},${job_id_7_6} --parsable 16_LIANA+_cell_multidonor.sh)
job_id_7_8=$(sbatch --dependency=afterok:${job_id_7_7} --parsable 16_LIANA+_cell_multidonor_extracell.sh)
job_id_7_9=$(sbatch --dependency=afterok:${job_id_7_8} --parsable 17_LIANA+_cell_type_spcific.sh)
job_id_7_10=$(sbatch --dependency=afterok:${job_id_7_9} --parsable 18_LIANA+_NMF.sh)
job_id_7_11=$(sbatch --dependency=afterok:${job_id_7_10} --parsable 19_GO_enrichment.sh)

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/

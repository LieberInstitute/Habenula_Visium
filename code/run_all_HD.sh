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
#   FICTURE with library-size-normalized inputs
################################################################################

cd $repo_dir/code/10_HD_bin_level/ficture_harmony

#-------------------------------------------------------------------------------
#   Main steps to run FICTURE
#-------------------------------------------------------------------------------

job_id_1_1=$(sbatch --parsable 01_build_spe.sh)
job_id_1_2=$(sbatch --dependency=afterok:${job_id_1_1} --parsable 02_normalized_input.sh)
job_id_1_3=$(sbatch --dependency=afterok:${job_id_1_2} --parsable 03_ficture_run.sh)
job_id_1_4=$(sbatch --dependency=afterok:${job_id_1_3} --parsable 03_ficture_plot.sh)
job_id_1_5=$(sbatch --dependency=afterok:${job_id_1_4} --parsable 04_spatula_join.sh)
job_id_1_6=$(sbatch --dependency=afterok:${job_id_1_5} --parsable 09_bin_level_merge_norm.sh)

#-------------------------------------------------------------------------------
#   Spatial registration
#-------------------------------------------------------------------------------

job_id_1_7=$(sbatch --dependency=afterok:${job_id_1_6} --parsable 10_registration_wrapper_norm.sh)
job_id_1_8=$(sbatch --dependency=afterok:${job_id_1_7} --parsable 11_cor_heatmap_norm.sh)

################################################################################
#   FICTURE with cleaningY batch-corrected inputs
################################################################################

#-------------------------------------------------------------------------------
#   Main steps to run FICTURE
#-------------------------------------------------------------------------------

job_id_2_1=$(sbatch --dependency=afterok:${job_id_1_1} --parsable 08_batcheffect_lm.sh)
job_id_2_2=$(sbatch --dependency=afterok:${job_id_2_1} --parsable 08_batcheffect_com.sh)
job_id_2_3=$(sbatch --dependency=afterok:${job_id_2_2} --parsable 08_reprepareinput_ficture.sh)
job_id_2_4=$(sbatch --dependency=afterok:${job_id_2_3} --parsable 08_rerun_ficture.sh)
job_id_2_5=$(sbatch --dependency=afterok:${job_id_2_4} --parsable 03_ficture_plot_cleany.sh)
job_id_2_6=$(sbatch --dependency=afterok:${job_id_2_5} --parsable 09_get_ficture_clusters.sh)
job_id_2_7=$(sbatch --dependency=afterok:${job_id_2_6} --parsable 09_bin_level_merge_cleany.sh)

#-------------------------------------------------------------------------------
#   Spatial registration
#-------------------------------------------------------------------------------

job_id_2_8=$(sbatch --dependency=afterok:${job_id_2_7} --parsable 10_registration_wrapper_cleany.sh)
job_id_2_9=$(sbatch --dependency=afterok:${job_id_2_8} --parsable 11_cor_heatmap_cleany.sh)

################################################################################
#   Create SpatialExperiments
################################################################################

#   Bin-level
cd $repo_dir/code/10_HD_bin_level
job_id_2_1=$(sbatch --parsable 01_build_spe.sh)
job_id_2_2=$(sbatch --dependency=afterok:${job_id_2_1} --parsable 02_QC.sh)

#   Cell-level
cd $repo_dir/code/09_HD_cell_level
job_id_2_3=$(sbatch --parsable 01_bin2cell.sh)
job_id_2_4=$(sbatch --dependency=afterok:${job_id_2_3} --parsable 02_build_spe_raw.sh)
job_id_2_5=$(sbatch --dependency=afterok:${job_id_2_4} --parsable 03_build_spe_QC.sh)
job_id_2_6=$(sbatch --dependency=afterok:${job_id_2_5} --parsable 04_HVG.sh)

################################################################################
#   Banksy
################################################################################

#   Finding SVGs
cd $repo_dir/code/10_HD_bin_level
job_id_3_1=$(sbatch --dependency=afterok:${job_id_2_6} --parsable 03_rasterize.sh)
job_id_3_2=$(sbatch --dependency=afterok:${job_id_3_1} --parsable 04_nnSVG.sh)
job_id_3_3=$(sbatch --dependency=afterok:${job_id_3_2} --parsable 05_gather_variable_genes.sh)

#   Running Banksy
cd $repo_dir/code/09_HD_cell_level
job_id_3_4=$(sbatch --dependency=afterok:${job_id_3_3} --parsable 05_banksy_embedding.sh)
job_id_3_5=$(sbatch --dependency=afterok:${job_id_3_4} --parsable 06_banksy_clustering.sh)

################################################################################
#   Spatial registration
################################################################################

#   Banksy
cd $repo_dir/code/09_HD_cell_level/registration_banksy
job_id_4_1=$(sbatch --dependency=afterok:${job_id_3_5} --parsable 01_registration_wrapper.sh)
job_id_4_2=$(sbatch --dependency=afterok:${job_id_4_1} --parsable 02_cor_heatmap.sh)

################################################################################
#   LIANA
################################################################################

cd $repo_dir/code/10_HD_bin_level/exploratory
job_id_5_1=$(sbatch --dependency=afterok:${job_id_4_2} --parsable 12_preprocess_anndata.sh)
job_id_5_2=$(sbatch --dependency=afterok:${job_id_5_1} --parsable 14_0_liana_prepare_inputs.sh)
job_id_5_3=$(sbatch --dependency=afterok:${job_id_5_2} --parsable 14_2_LIANA+_preprocess.sh)
job_id_5_4=$(sbatch --dependency=afterok:${job_id_5_3} --parsable 14_3_LIANA+_preprocss_extracellular.sh)
job_id_5_5=$(sbatch --dependency=afterok:${job_id_5_4} --parsable 15_1_LIANA+_cell.sh)
job_id_5_6=$(sbatch --dependency=afterok:${job_id_5_4} --parsable 15_1_LIANA+_extracell.sh)
job_id_5_7=$(sbatch --dependency=afterok:${job_id_5_5},${job_id_5_6} --parsable 16_LIANA+_cell_multidonor.sh)
job_id_5_8=$(sbatch --dependency=afterok:${job_id_5_7} --parsable 16_LIANA+_cell_multidonor_extracell.sh)
job_id_5_9=$(sbatch --dependency=afterok:${job_id_5_8} --parsable 17_LIANA+_cell_type_spcific.sh)
job_id_5_10=$(sbatch --dependency=afterok:${job_id_5_9} --parsable 18_LIANA+_NMF.sh)
job_id_5_11=$(sbatch --dependency=afterok:${job_id_5_10} --parsable 19_GO_enrichment.sh)

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/

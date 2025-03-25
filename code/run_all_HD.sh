#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=run_all
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH -o ../run_all_hd.txt
#SBATCH -e ../run_all_hd.txt

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

################################################################################
#   FICTURE with library-size-normalized inputs
################################################################################

cd 10_HD_bin_level/ficture_harmony
job_id_1_1=$(sbatch --parsable 01_build_spe.sh)
job_id_2_1=$(sbatch --dependency=afterok:${job_id_1_1} --parsable 02_normalized_input.sh)
job_id_3_1=$(sbatch --dependency=afterok:${job_id_2_1} --parsable 03_ficture_run.sh)
job_id_4_1=$(sbatch --dependency=afterok:${job_id_3_1} --parsable 04_spatula_join.sh)
job_id_5_1=$(sbatch --dependency=afterok:${job_id_4_1} --parsable 05_bin2cell.sh)
job_id_6_1=$(sbatch --dependency=afterok:${job_id_5_1} --parsable 06_merge_and_plot.sh)
cd ../..

################################################################################
#   Create SpatialExperiments
################################################################################

job_id_1=$(sbatch --parsable 01_build_spe.sh)
job_id_2=$(sbatch --dependency=afterok:${job_id_1} --parsable 02_QC.sh)
cd 09_HD_cell_level
job_id_3=$(sbatch --parsable 01_bin2cell.sh)
job_id_4=$(sbatch --dependency=afterok:${job_id_3} --parsable 02_build_spe_raw.sh)
job_id_5=$(sbatch --dependency=afterok:${job_id_4} --parsable 03_build_spe_QC.sh)
cd ..

################################################################################
#   Banksy
################################################################################

#   Finding SVGs
cd 10_HD_bin_level
job_id_6=$(sbatch --dependency=afterok:${job_id_5} --parsable 06_rasterize.sh)
job_id_7=$(sbatch --dependency=afterok:${job_id_6} --parsable 07_nnSVG.sh)
job_id_8=$(sbatch --dependency=afterok:${job_id_7} --parsable 08_gather_variable_genes.sh)

#   Running Banksy
cd ../09_HD_cell_level
job_id_9=$(sbatch --dependency=afterok:${job_id_8} --parsable 08_banksy_embedding.sh)
job_id_10=$(sbatch --dependency=afterok:${job_id_9} --parsable 10_banksy_clustering.sh)

################################################################################
#   Spatial registration
################################################################################

#   Banksy
cd registration_banksy
job_id_11=$(sbatch --dependency=afterok:${job_id_10} --parsable 01_registration_wrapper.sh)
job_id_12=$(sbatch --dependency=afterok:${job_id_11} --parsable 02_cor_heatmap.sh)
cd ../../10_HD_bin_level/ficture_harmony

#   FICTURE
job_id_13=$(sbatch --dependency=afterok:${job_id_6_1} --parsable 07_registration_wrapper.sh)

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/

#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=30G
#SBATCH --job-name=06_combine_ficture_results
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../../processed-data/10_HD_bin_level/logs/all_sample_002um/06_combine_ficture_results_%a.log
#SBATCH -e ../../../processed-data/10_HD_bin_level/logs/all_sample_002um/06_combine_ficture_results_%a.log
#SBATCH --array=2-25%10

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

ml spatula/f0e9936

repo_dir=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium
sample_id_path=$repo_dir/raw-data/sample_info/hd_sample_list.txt

#rerun the join-pixel-tsv
out_dir=$repo_dir/processed-data/10_HD_bin_level/ficture/outputs/all_samples/k_$SLURM_ARRAY_TASK_ID/analysis/nF${SLURM_ARRAY_TASK_ID}.d_${SLURM_ARRAY_TASK_ID}
in_tsv=$repo_dir/processed-data/10_HD_bin_level/ficture/inputs/all_samples/transcripts_moved_with_barcodes_sorted.tsv.gz

#   Sort FICTURE output by major axis
(gzip -cd $out_dir/nF${SLURM_ARRAY_TASK_ID}.d_${SLURM_ARRAY_TASK_ID}.decode.prj_${SLURM_ARRAY_TASK_ID}.r_4_5.pixel.sorted.tsv.gz \
    | head | grep ^#; \
    gzip -cd $out_dir/nF${SLURM_ARRAY_TASK_ID}.d_${SLURM_ARRAY_TASK_ID}.decode.prj_${SLURM_ARRAY_TASK_ID}.r_4_5.pixel.sorted.tsv.gz \
    | grep -v ^# | sort -S 1G -gk2) \
    | gzip -c > $out_dir/nF${SLURM_ARRAY_TASK_ID}.d_${SLURM_ARRAY_TASK_ID}.decode.prj_${SLURM_ARRAY_TASK_ID}.r_4_5.pixel.sorted_by_major_axis.tsv.gz

#   Write transcript-level output
spatula join-pixel-tsv \
    --mol-tsv $in_tsv \
    --pix-prefix-tsv factor_,$out_dir/nF${SLURM_ARRAY_TASK_ID}.d_${SLURM_ARRAY_TASK_ID}.decode.prj_${SLURM_ARRAY_TASK_ID}.r_4_5.pixel.sorted_by_major_axis.tsv.gz \
    --out-prefix $out_dir/transcripts_ficture_joined_moved_with_barcodes \
    --out-max-k 3 \
    --out-max-p 3 \
    --mu-scale 1

echo "**** Job ends ****"
date

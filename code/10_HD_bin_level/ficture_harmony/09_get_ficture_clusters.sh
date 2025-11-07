#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=30G
#SBATCH --job-name=09_get_ficture_clusters
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../../processed-data/10_HD_bin_level/new_samples/ficture_harmony/logs/09_get_ficture_clusters_%a.txt
#SBATCH -e ../../../processed-data/10_HD_bin_level/new_samples/ficture_harmony/logs/09_get_ficture_clusters_%a.txt
#SBATCH --array=3-29%10

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

ml spatula/f0e9936

a=$SLURM_ARRAY_TASK_ID
#rerun the join-pixel-tsv
repo_dir=$(git rev-parse --show-toplevel)
out_dir=$repo_dir/processed-data/10_HD_bin_level/new_samples/ficture_harmony/ficture_outputs/cleany/k_${a}/analysis/nF${a}.d_12
in_tsv=$repo_dir/processed-data/10_HD_bin_level/new_samples/ficture_harmony/ficture_inputs/cleany/input.tsv.gz

#   Sort FICTURE output by major axis
(gzip -cd $out_dir/nF${a}.d_12.decode.prj_12.r_4_5.pixel.sorted.tsv.gz \
    | head | grep ^#; \
    gzip -cd $out_dir/nF${a}.d_12.decode.prj_12.r_4_5.pixel.sorted.tsv.gz \
    | grep -v ^# | sort -S 1G -gk2) \
    | gzip -c > $out_dir/nF${a}.d_12.decode.prj_12.r_4_5.pixel.sorted_by_major_axis.tsv.gz

#   Write transcript-level output
spatula join-pixel-tsv \
    --mol-tsv $in_tsv \
    --pix-prefix-tsv factor_,$out_dir/nF${a}.d_12.decode.prj_12.r_4_5.pixel.sorted_by_major_axis.tsv.gz \
    --out-prefix $out_dir/transcripts_ficture_joined \
    --out-max-k 3 \
    --out-max-p 3 \
    --mu-scale 1

echo "**** Job ends ****"
date

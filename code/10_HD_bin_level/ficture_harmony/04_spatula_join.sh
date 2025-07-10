#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=30G
#SBATCH --job-name=04_spatula_join
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../../processed-data/10_HD_bin_level/probe_fix/ficture_harmony/logs/04_spatula_join_%a.txt
#SBATCH -e ../../../processed-data/10_HD_bin_level/probe_fix/ficture_harmony/logs/04_spatula_join_%a.txt
#SBATCH --array=35,37,39%10

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

ml spatula/f0e9936

repo_dir=$(git rev-parse --show-toplevel)

#   Path definitions
out_dir=$repo_dir/processed-data/10_HD_bin_level/probe_fix/ficture_harmony/ficture_outputs/normalized/k_$SLURM_ARRAY_TASK_ID/analysis/nF${SLURM_ARRAY_TASK_ID}.d_12
in_tsv=$repo_dir/processed-data/10_HD_bin_level/probe_fix/ficture_harmony/ficture_inputs/normalized/input.tsv.gz

#   Sort FICTURE output by major axis
(gzip -cd $out_dir/nF${SLURM_ARRAY_TASK_ID}.d_12.decode.prj_12.r_4_5.pixel.sorted.tsv.gz \
    | head | grep ^#; \
    gzip -cd $out_dir/nF${SLURM_ARRAY_TASK_ID}.d_12.decode.prj_12.r_4_5.pixel.sorted.tsv.gz \
    | grep -v ^# | sort -S 1G -gk2) \
    | gzip -c > $out_dir/nF${SLURM_ARRAY_TASK_ID}.d_12.decode.prj_12.r_4_5.pixel.sorted_by_major_axis.tsv.gz

#   Write transcript-level output
if [[ $SLURM_ARRAY_TASK_ID -eq 2 ]]; then
    spatula join-pixel-tsv \
        --mol-tsv $in_tsv \
        --pix-prefix-tsv factor_,$out_dir/nF${SLURM_ARRAY_TASK_ID}.d_12.decode.prj_12.r_4_5.pixel.sorted_by_major_axis.tsv.gz \
        --out-prefix $out_dir/normalized_joined_input \
        --out-max-k 2 \
        --out-max-p 2 \
        --mu-scale 1
else
    spatula join-pixel-tsv \
        --mol-tsv $in_tsv \
        --pix-prefix-tsv factor_,$out_dir/nF${SLURM_ARRAY_TASK_ID}.d_12.decode.prj_12.r_4_5.pixel.sorted_by_major_axis.tsv.gz \
        --out-prefix $out_dir/normalized_joined_input \
        --out-max-k 3 \
        --out-max-p 3 \
        --mu-scale 1
fi

echo "**** Job ends ****"
date

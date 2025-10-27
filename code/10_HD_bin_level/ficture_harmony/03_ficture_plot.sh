#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=5G
#SBATCH --job-name=03_ficture_plot
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH -o ../../../processed-data/10_HD_bin_level/new_samples/ficture_harmony/logs/03_ficture_plot_%a.txt
#SBATCH -e ../../../processed-data/10_HD_bin_level/new_samples/ficture_harmony/logs/03_ficture_plot_%a.txt
#SBATCH --array=4-40:2,70,100%15

#   Default plots produced by 'ficture run_together' have too much empty (black)
#   space. Reproduce these plots with visually preferable settings (for
#   library-size normalized results)

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load visium_hd/1.0

repo_dir=$(git rev-parse --show-toplevel)
out_dir=$repo_dir/processed-data/10_HD_bin_level/new_samples/ficture_harmony/ficture_outputs/normalized/k_${SLURM_ARRAY_TASK_ID}
plot_dir=$repo_dir/plots/10_HD_bin_level/new_samples/ficture_harmony/normalized

mkdir -p $plot_dir

ficture plot_pixel_full \
    --input $out_dir/analysis/nF${SLURM_ARRAY_TASK_ID}.d_12/nF${SLURM_ARRAY_TASK_ID}.d_12.decode.prj_12.r_4_5.pixel.sorted.tsv.gz \
    --color_table $out_dir/analysis/nF${SLURM_ARRAY_TASK_ID}.d_12/figure/nF${SLURM_ARRAY_TASK_ID}.d_12.rgb.tsv \
    --output $plot_dir/k_${SLURM_ARRAY_TASK_ID}.png \
    --plot_um_per_pixel 1.5 \
    --full

echo "**** Job ends ****"
date

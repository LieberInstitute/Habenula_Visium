#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=04_ficture_plot
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../processed-data/10_HD_bin_level/logs/04_ficture_plot_all_samples.log
#SBATCH -e ../../processed-data/10_HD_bin_level/logs/04_ficture_plot_all_samples.log

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"

module load visium_hd/1.0

out_dir=$repo_dir/processed-data/10_HD_bin_level/ficture/outputs/all_samples
plot_dir=$repo_dir/plots/10_HD_bin_level/ficture/all_sample_002um

mkdir -p $plot_dir

ficture plot_pixel_full \
    --input $out_dir/analysis/nF12.d_12/nF12.d_12.decode.prj_12.r_4_5.pixel.sorted.tsv.gz \
    --color_table $out_dir/analysis/nF12.d_12/figure/nF12.d_12.rgb.tsv \
    --output $plot_dir/clean_clusters_all_samples.png \
    --plot_um_per_pixel 1.5 \
    --full

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.2
## available from http://research.libd.org/slurmjobs/

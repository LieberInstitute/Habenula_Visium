#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=07_ficture_plot
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH -o ../../../processed-data/10_HD_bin_level/logs/all_sample_002um/07_ficture_plot_%a.log
#SBATCH -e ../../../processed-data/10_HD_bin_level/logs/all_sample_002um/07_ficture_plot_%a.log
#SBATCH --array=1-5%5

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"

ml ficture/dev_a455e5c
ml spatula/f0e9936

repo_dir=$(git rev-parse --show-toplevel)
sample_id_path=$repo_dir/raw-data/sample_info/hd_sample_list.txt
this_sample=$(awk "NR==${SLURM_ARRAY_TASK_ID}" $sample_id_path)

temp_dir=$MYSCRATCH
out_dir=$repo_dir/processed-data/10_HD_bin_level/ficture/outputs/all_samples
plot_dir=$repo_dir/plots/10_HD_bin_level/ficture

mkdir -p $plot_dir

#   Subset the output to just this sample for plotting
zcat $out_dir/analysis/nF12.d_12/nF12.d_12.decode.prj_12.r_4_5.pixel.sorted.tsv.gz \
    | sed -E 's|SIZE_X=([0-9]+)|SIZE_X=7000|' \
    | awk -F '\t' -v low_cutoff=$((7000 * ($SLURM_ARRAY_TASK_ID - 1))) '
        BEGIN {OFS="\t"}
        NR<=4 {print $0}
        NR>4 && ($2 / 100 >= low_cutoff) && ($2 / 100 < low_cutoff + 7000) {
            $2 = $2 - low_cutoff * 100; print $0
        }' \
    | gzip -c > $temp_dir/${this_sample}.tsv.gz

ficture plot_pixel_full \
    --input $temp_dir/${this_sample}.tsv.gz \
    --color_table $out_dir/analysis/nF12.d_12/figure/nF12.d_12.rgb.tsv \
    --output $plot_dir/clean_clusters_${this_sample}.png \
    --plot_um_per_pixel 1.5 \
    --full

rm $temp_dir/${this_sample}.tsv.gz

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.2
## available from http://research.libd.org/slurmjobs/

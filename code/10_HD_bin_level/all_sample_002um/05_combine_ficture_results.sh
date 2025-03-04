#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=50G
#SBATCH --job-name=05_combine_ficture_results
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../../processed-data/10_HD_bin_level/logs/all_sample_002um/05_combine_ficture_results.log
#SBATCH -e ../../../processed-data/10_HD_bin_level/logs/all_sample_002um/05_combine_ficture_results.log

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

ml ficture/dev_a455e5c
ml spatula/f0e9936

repo_dir=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium
sample_id_path=$repo_dir/raw-data/sample_info/hd_sample_list.txt
temp_dir=$MYSCRATCH

data_dir=$repo_dir/processed-data/01_spaceranger/$this_sample/outs/binned_outputs/square_002um
out_dir=$repo_dir/processed-data/10_HD_bin_level/ficture/inputs/$this_sample

# Merge 5 sample together
merged_output=$repo_dir/processed-data/10_HD_bin_level/ficture/inputs/all_samples/transcripts_moved_with_barcodes_sorted.tsv.gz

for i in $(seq 1 5); do
    this_sample=$(awk "NR==$i" "$sample_id_path")
    out_dir="$repo_dir/processed-data/10_HD_bin_level/ficture/inputs/$this_sample"

    if [ "$i" -eq 1 ]; then
        zcat "$out_dir/transcripts_moved_with_barcodes_sorted.tsv.gz" > "temp_file"
    else
        zcat "$out_dir/transcripts_moved_with_barcodes_sorted.tsv.gz" | tail -n +2 >> "temp_file"
    fi
done

gzip -c "temp_file" > "$merged_output"


#rerun the join-pixel-tsv
out_dir=$repo_dir/processed-data/10_HD_bin_level/ficture/outputs/all_samples/analysis/nF12.d_12
in_tsv=$repo_dir/processed-data/10_HD_bin_level/ficture/inputs/all_samples/transcripts_moved_with_barcodes_sorted.tsv.gz

#   Sort FICTURE output by major axis
(gzip -cd $out_dir/nF12.d_12.decode.prj_12.r_4_5.pixel.sorted.tsv.gz \
    | head | grep ^#; \
    gzip -cd $out_dir/nF12.d_12.decode.prj_12.r_4_5.pixel.sorted.tsv.gz \
    | grep -v ^# | sort -S 1G -gk2) \
    | gzip -c > $out_dir/nF12.d_12.decode.prj_12.r_4_5.pixel.sorted_by_major_axis.tsv.gz

#   Write transcript-level output
spatula join-pixel-tsv \
    --mol-tsv $in_tsv \
    --pix-prefix-tsv factor_,$out_dir/nF12.d_12.decode.prj_12.r_4_5.pixel.sorted_by_major_axis.tsv.gz \
    --out-prefix $out_dir/transcripts_ficture_joined_moved_with_barcodes \
    --max-dist-um 2 \
    --bin-um 4 \
    --out-max-k 3 \
    --out-max-p 3 \
    --mu-scale 1

echo "**** Job ends ****"
date
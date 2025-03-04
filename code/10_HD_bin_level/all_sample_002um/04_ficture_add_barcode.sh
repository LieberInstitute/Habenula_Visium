#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=04_ficture_add_barcode.sh
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../../processed-data/10_HD_bin_level/logs/all_sample_002um/04_ficture_add_barcode_%a.log
#SBATCH -e ../../../processed-data/10_HD_bin_level/logs/all_sample_002um/04_ficture_add_barcode_%a.log
#SBATCH --array=1-5%5

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
this_sample=$(awk "NR==${SLURM_ARRAY_TASK_ID}" $sample_id_path)

data_dir=$repo_dir/processed-data/01_spaceranger/$this_sample/outs/binned_outputs/square_002um
out_dir=$repo_dir/processed-data/10_HD_bin_level/ficture/inputs/$this_sample

# Set file paths
matrix="$data_dir/raw_feature_bc_matrix/matrix.mtx.gz"
barcodes="$data_dir/raw_feature_bc_matrix/barcodes.tsv.gz"
transcripts="$out_dir/transcripts.unsorted.tsv.gz"
output="$out_dir/transcripts_with_barcodes.tsv.gz"

# 1. Read barcodes and store them in an indexed file
zcat "$barcodes" | awk '{print NR, $0}' > $out_dir/barcodes_index.txt

# 2. Read matrix.mtx.gz, match gene-cell pairs with barcodes
zcat "$matrix" | tail -n +4 | awk 'NR==FNR {barcode[$1] = $2; next} {print $0, barcode[$2]}' $out_dir/barcodes_index.txt - > $out_dir/matrix_with_barcodes.tsv

# 3. Merge barcodes into transcripts.unsorted.tsv.gz
(echo -e "$(zcat "$transcripts" | head -1)\tbarcode"
 paste <(zcat "$transcripts" | tail -n +2) <(cut -d' ' -f4 $out_dir/matrix_with_barcodes.tsv)
) | gzip > "$output"

zcat $output | head -10
zcat $output | tail -10

# Remove temporary files
rm $out_dir/barcodes_index.txt $out_dir/matrix_with_barcodes.tsv

echo "✅ Merging completed! Output file: $output"

output="$out_dir/transcripts_with_barcodes.tsv.gz"

# Sort this transcript
gzip -cd "$output" | (
    head -1
    tail -n +2 | sort -S 1G -g -k1,1
) | gzip -c > "$out_dir/transcripts_with_barcodes_sorted.tsv.gz"

zcat $out_dir/transcripts_with_barcodes_sorted.tsv.gz | head -10
zcat $out_dir/transcripts_with_barcodes_sorted.tsv.gz | tail -10

# Move the 5 samples
xmin=$(grep "xmin" $out_dir/minmax.tsv | cut -f 2)
ymin=$(grep "ymin" $out_dir/minmax.tsv | cut -f 2)
zcat $out_dir/transcripts_with_barcodes_sorted.tsv.gz | \
    awk \
    -v min_x="$xmin" \
    -v min_y="$ymin" \
    -v offset="$((7000 * (${SLURM_ARRAY_TASK_ID} - 1)))" \
    'NR==1 {print; next} { $1 = $1 + offset - min_x; $2 = $2 - min_y; print }' OFS="\t" \
    | gzip > $out_dir/transcripts_moved_with_barcodes_sorted.tsv.gz

echo "✅ Moving completed! Output file: $out_dir/transcripts_moved_with_barcodes_sorted.tsv.gz"

echo "**** Job ends ****"
date
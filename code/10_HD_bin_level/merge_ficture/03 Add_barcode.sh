#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=50G
#SBATCH --job-name=04_ficture_03_Add_barcode.sh
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../processed-data/10_HD_bin_level/logs/04_ficture_03_Add_barcode_%a.log
#SBATCH -e ../../processed-data/10_HD_bin_level/logs/04_ficture_03_Add_barcode_%a.log
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

data_dir=$repo_dir/processed-data/01_spaceranger/$this_sample/outs/binned_outputs/square_008um
out_dir=$repo_dir/processed-data/10_HD_bin_level/ficture/inputs/square_008um/$this_sample

# Set file paths
matrix="$data_dir/raw_feature_bc_matrix/matrix.mtx.gz"
barcodes="$data_dir/raw_feature_bc_matrix/barcodes.tsv.gz"
transcripts="$out_dir/transcripts.unsorted.tsv.gz"
output="$out_dir/transcripts_with_barcodes.tsv.gz"

# 1. Read barcodes and store them in an indexed file
zcat "$barcodes" | awk '{print NR, $0}' > barcodes_index.txt

# 2. Read matrix.mtx.gz, match gene-cell pairs with barcodes
zcat "$matrix" | tail -n +4 | awk 'NR==FNR {barcode[$1] = $2; next} {print $0, barcode[$2]}' barcodes_index.txt - > matrix_with_barcodes.tsv

# 3. Merge barcodes into transcripts.unsorted.tsv.gz
(echo -e "$(zcat "$transcripts" | head -1)\tbarcode"
 paste <(zcat "$transcripts" | tail -n +2) <(cut -d' ' -f4 matrix_with_barcodes.tsv)
) | gzip > "$output"

zcat $output | head -10

# Remove temporary files
rm barcodes_index.txt matrix_with_barcodes.tsv

echo "✅ Merging completed! Output file: $output"

# Sort this transcript
(gzip -cd $output \
    | head -1; gzip -cd $output \
    | tail -n +2 | sort -S 1G -gk1) \
    | gzip -c > $out_dir/transcripts_with_barcodes_sorted.tsv.gz

#rerun the join-pixel-tsv
out_dir=$repo_dir/processed-data/10_HD_bin_level/ficture/outputs/$this_sample/analysis/nF12.d_12
in_tsv=$repo_dir/processed-data/10_HD_bin_level/ficture/inputs/square_008um/$this_sample/transcripts_with_barcodes_sorted.tsv.gz

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
    --out-prefix $out_dir/transcripts_ficture_joined \
    --max-dist-um 2 \
    --bin-um 4 \
    --out-max-k 3 \
    --out-max-p 3 \
    --mu-scale 1

echo "**** Job ends ****"
date
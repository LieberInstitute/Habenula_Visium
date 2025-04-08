#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=05_combine_barcodes
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../../processed-data/10_HD_bin_level/logs/all_sample_002um/05_combine_barcodes.log
#SBATCH -e ../../../processed-data/10_HD_bin_level/logs/all_sample_002um/05_combine_barcodes.log

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

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

    #   Append each TSV to an ongoing merged file. Add sample ID as a constant
    #   column
    if [ "$i" -eq 1 ]; then
        zcat "$out_dir/transcripts_moved_with_barcodes_sorted.tsv.gz" \
            | awk -F '\t' -v this_sample=$this_sample 'NR==1 {print $0 "\tsample_id"} NR>1 {print $0 "\t" this_sample}' \
            > $temp_dir/temp.tsv
    else
        zcat "$out_dir/transcripts_moved_with_barcodes_sorted.tsv.gz" \
            | awk -F '\t' -v this_sample=$this_sample '{print $0 "\t" this_sample}' \
            | tail -n +2 \
            >> $temp_dir/temp.tsv
    fi
done

#   Compress merged file
gzip -c $temp_dir/temp.tsv > "$merged_output"
rm $temp_dir/temp.tsv

echo "**** Job ends ****"
date
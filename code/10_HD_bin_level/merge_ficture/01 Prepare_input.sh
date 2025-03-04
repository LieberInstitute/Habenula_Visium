#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=50G
#SBATCH --job-name=04_ficture_01_Prepare_input
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../processed-data/10_HD_bin_level/logs/04_ficture_01_Prepare_input_%a.log
#SBATCH -e ../../processed-data/10_HD_bin_level/logs/04_ficture_01_Prepare_input_%a.log
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

echo "Processing sample: $this_sample"

data_dir=$repo_dir/processed-data/01_spaceranger/$this_sample/outs/binned_outputs/square_008um
out_dir=$repo_dir/processed-data/10_HD_bin_level/ficture/inputs/square_008um/$this_sample

mkdir -p $out_dir

#   Get spatial coordinates as a CSV
echo "Converting spatial coords to CSV..."
parquet-tools csv $data_dir/spatial/tissue_positions.parquet \
    | gzip -c > $temp_dir/${this_sample}_tissue_positions.csv.gz

microns_per_pixel=$(
    grep microns_per_pixel $data_dir/spatial/scalefactors_json.json \
        | cut -d ":" -f 2 \
        | tr -d ", "
)

spatula convert-sge \
    --in-sge $data_dir/raw_feature_bc_matrix \
    --pos $temp_dir/${this_sample}_tissue_positions.csv.gz \
    --units-per-um $(python -c "print(1/${microns_per_pixel})") \
    --colnames-count Count \
    --out-tsv $out_dir \
    --icols-mtx 1

## Sort the unsorted output file by the X-coordinate
(gzip -cd $out_dir/transcripts.unsorted.tsv.gz \
    | head -1; gzip -cd $out_dir/transcripts.unsorted.tsv.gz \
    | tail -n +2 | sort -S 1G -gk1) \
    | gzip -c > $out_dir/transcripts.sorted.tsv.gz

echo "**** Job ends ****"
date
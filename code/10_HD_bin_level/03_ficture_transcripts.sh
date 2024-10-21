#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=5G
#SBATCH --job-name=03_ficture_transcripts
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../processed-data/10_HD_bin_level/logs/03_ficture_transcripts.log
#SBATCH -e ../../processed-data/10_HD_bin_level/logs/03_ficture_transcripts.log

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load ficture/0.0.3.1

repo_dir=$(git rev-parse --show-toplevel)
temp_dir=$MYSCRATCH
this_sample=H1-W369TJK_D1_9090

data_dir=$repo_dir/processed-data/01_spaceranger/$this_sample/outs/binned_outputs/square_002um
out_dir=$repo_dir/processed-data/10_HD_bin_level/ficture/inputs/$this_sample

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

## This script was made using slurmjobs version 1.2.2
## available from http://research.libd.org/slurmjobs/

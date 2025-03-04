#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=03_ficture_run
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../../processed-data/10_HD_bin_level/logs/all_sample_002um/03_ficture_run.log
#SBATCH -e ../../../processed-data/10_HD_bin_level/logs/all_sample_002um/03_ficture_run.log

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

repo_dir=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium
sample_id_path=$repo_dir/raw-data/sample_info/hd_sample_list.txt

mkdir -p $repo_dir/processed-data/10_HD_bin_level/ficture/inputs/all_samples

#  Transcripts.moved.sorted.tsv.gz
merged_output=$repo_dir/processed-data/10_HD_bin_level/ficture/inputs/all_samples/transcripts.merged.sorted.tsv.gz

for i in $(seq 1 5); do
    this_sample=$(awk "NR==$i" "$sample_id_path")
    out_dir="$repo_dir/processed-data/10_HD_bin_level/ficture/inputs/$this_sample"

    if [ "$i" -eq 1 ]; then
        zcat "$out_dir/transcripts.moved.sorted.tsv.gz" > "temp_file"
    else
        zcat "$out_dir/transcripts.moved.sorted.tsv.gz" | tail -n +2 >> "temp_file"
    fi
done

gzip -c "temp_file" > "$merged_output"

# minmax.tsv
awk 'NR>1 {
    if($1+0==$1 && $2+0==$2) {
        if(NR==2 || $1<xmin) xmin=$1;
        if(NR==2 || $1>xmax) xmax=$1;
        if(NR==2 || $2<ymin) ymin=$2;
        if(NR==2 || $2>ymax) ymax=$2;
    }
} 
END {
    print "xmin\t" xmin "\nxmax\t" xmax "\nymin\t" ymin "\nymax\t" ymax;
}' temp_file > $repo_dir/processed-data/10_HD_bin_level/ficture/inputs/all_samples/merged_minmax.tsv

rm "temp_file"  


#   Path definitions
in_dir=$repo_dir/processed-data/10_HD_bin_level/ficture/inputs/all_samples
out_dir=$repo_dir/processed-data/10_HD_bin_level/ficture/outputs/all_samples

mkdir -p $out_dir

#   Run full FICTURE pipeline
ficture run_together \
    --in-tsv $in_dir/transcripts.merged.sorted.tsv.gz \
    --in-minmax $in_dir/merged_minmax.tsv \
    --out-dir $out_dir \
    --mu-scale 1 \
    --major-axis X \
    --all

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.2
## available from http://research.libd.org/slurmjobs/

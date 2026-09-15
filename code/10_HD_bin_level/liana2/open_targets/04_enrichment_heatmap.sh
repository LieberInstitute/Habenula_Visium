#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=04_enrichment_heatmap.py
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../../../processed-data/10_HD_bin_level/no_secondary/liana2/open_targets/logs/04_enrichment_heatmap.txt
#SBATCH -e ../../../../processed-data/10_HD_bin_level/no_secondary/liana2/open_targets/logs/04_enrichment_heatmap.txt

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load the R module
module load liana_plus/1.7.1

## List current modules for reproducibility
module list

bandwidth="5000"
disease_name="Substance_dependence"

# python3 '04_enrichment_heatmap.py' --bandwidth "$bandwidth" --disease "$disease_name"
python3 '04_enrichment_heatmap_either.py' --bandwidth "$bandwidth" --disease "$disease_name" --min_donors 4

bandwidth="5000"
disease_name="MDD"

# python3 '04_enrichment_heatmap.py' --bandwidth "$bandwidth" --disease "$disease_name"
python3 '04_enrichment_heatmap_either.py' --bandwidth "$bandwidth" --disease "$disease_name" --min_donors 4

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/


#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=06_geneset_enrichment.py
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../../../processed-data/10_HD_bin_level/no_secondary/liana2/open_targets/logs/06_geneset_enrichment.txt
#SBATCH -e ../../../../processed-data/10_HD_bin_level/no_secondary/liana2/open_targets/logs/06_geneset_enrichment.txt

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load the Python module
module load liana_plus/1.7.1

## List current modules for reproducibility
module list

SCRIPT_DIR="/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/code/10_HD_bin_level/liana2/open_targets"

python3 "${SCRIPT_DIR}/06_geneset_enrichment.py" \
    --gene_sets "${SCRIPT_DIR}/trios_intersect_gene_sets.csv" \
    --disease MDD Substance_dependence \
    --background "${SCRIPT_DIR}/background_genes_protein_coding.csv"

echo "**** Job ends ****"
date


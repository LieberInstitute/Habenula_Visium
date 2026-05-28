#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=05_MAGMA_gene_set
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../../../processed-data/10_HD_bin_level/no_secondary/cell_environment/MAGMA/logs/05_MAGMA_gene_set_%a.txt
#SBATCH -e ../../../../processed-data/10_HD_bin_level/no_secondary/cell_environment/MAGMA/logs/05_MAGMA_gene_set_%a.txt
#SBATCH --array=1-16%16

#   Just the gene-set analysis step of MAGMA

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load magma/1.10
module list

all_gwas=(MDD panic SCZ SUD2020 AUD CUD ext_cannabis lifetime_cannabis OUD SUD2 compulsive internalizing neurodev p_factor SCZ_BPD SUD3)
gwas=${all_gwas[$(($SLURM_ARRAY_TASK_ID - 1))]}

repo_dir=$(git rev-parse --show-toplevel)
out_dir=${repo_dir}/processed-data/10_HD_bin_level/no_secondary/cell_environment/MAGMA/$gwas
gene_set_dir=${repo_dir}/processed-data/10_HD_bin_level/no_secondary/cell_environment/MAGMA/gene_sets

echo "Processing GWAS ${gwas}"

#   Gene set analysis step
for k in $(seq 3 10) 20; do
    mkdir -p $out_dir/$gwas/k$k
    magma \
        --gene-results $out_dir/$gwas.genes.raw \
        --set-annot $gene_set_dir/k$k.tsv gene-col=gene_id set-col=set_id \
        --out $out_dir/$gwas/k$k
done

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/

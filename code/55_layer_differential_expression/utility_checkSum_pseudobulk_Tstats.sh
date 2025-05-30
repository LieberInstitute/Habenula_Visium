#!/bin/bash
#SBATCH --partition=katun
#SBATCH --job-name=checksum_comparison
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH --mem=2GB
# SBATCH -o /dev/null
# SBATCH -e /dev/null
# SBATCH --mail-type=ALL

# Loop over the range from 2 to 28

echo "==========================="
echo "      TEST 1               "
echo "==========================="

echo "Checksum Comparison Results: 05_layer_differential_expression" > checksum_comparison_pseudobulk.txt
echo "===========================" >> checksum_comparison_pseudobulk.txt
echo "Test2: Comparing pseudobulk rds object" >> checksum_comparison_pseudobulk.txt

# rm -f slurm-*.out

dir1="../../processed-data/05_layer_differential_expression"
dir2="../../processed-data/05_layer_differential_expression/old"

# Loop over the range from 2 to 28
for i in $(seq -f "%02g" 2 28); do
  # Add a separator for each comparison
  echo "Comparing k${i}:" >> checksum_comparison_pseudobulk.txt
  sha256sum ${dir1}/sce_pseudo_BayesSpace_k${i}.rds ${dir2}/sce_pseudo_BayesSpace_k${i}.rds >> checksum_comparison_pseudobulk.txt
  echo "--------------------------------" >> checksum_comparison_pseudobulk.txt
done

# dumy example to illustrate not identical object
echo "dumy example to illustrate not identical object k02 vs K03:" >> checksum_comparison_pseudobulk.txt
sha256sum ${dir1}/sce_pseudo_BayesSpace_k02.rds ${dir2}/sce_pseudo_BayesSpace_k03.rds >> checksum_comparison_pseudobulk.txt
echo "--------------------------------" >> checksum_comparison_pseudobulk.txt


echo "==========================="
echo "      TEST 2               "
echo "==========================="

echo "Now I am comparing the enrichment object"

echo "Checksum Comparison Results: 05_layer_differential_expression" > checksum_comparison_enrichment.txt
echo "===========================" >> checksum_comparison_enrichment.txt
echo "Test2: Comparing enrichment Rdata object" >> checksum_comparison_enrichment.txt

# rm -f slurm-*.out

dir1="../../processed-data/05_layer_differential_expression/modeling_results_BS"
dir2="../../processed-data/05_layer_differential_expression/modeling_results_BS_old"

# Loop over the range from 2 to 28
for i in $(seq -f "%02g" 2 28); do
  # Add a separator for each comparison
  echo "Comparing enrichment k${i}:" >> checksum_comparison_enrichment.txt
  sha256sum ${dir1}/modeling_results_BayesSpace_k${i}.Rdata ${dir2}/modeling_results_BayesSpace_k${i}.Rdata >> checksum_comparison_enrichment.txt
  echo "--------------------------------" >> checksum_comparison_enrichment.txt
done

echo "Done!"


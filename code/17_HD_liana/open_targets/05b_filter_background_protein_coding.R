#!/usr/bin/env Rscript

# Filter background_genes.csv to protein-coding genes using biomaRt.
# Outputs: background_genes_protein_coding.csv

library(biomaRt)

script_dir <- "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/code/10_HD_bin_level/liana2/open_targets"

# Load background genes
bg <- read.csv(file.path(script_dir, "background_genes.csv"))
cat("Total background genes:", nrow(bg), "\n")

# Query biomaRt for protein-coding genes
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

pc_genes <- getBM(
    attributes = c("hgnc_symbol", "gene_biotype"),
    filters = "biotype",
    values = "protein_coding",
    mart = ensembl
)

pc_symbols <- unique(pc_genes$hgnc_symbol[pc_genes$hgnc_symbol != ""])
cat("Protein-coding gene symbols from biomaRt:", length(pc_symbols), "\n")

# Filter background to protein-coding
bg_pc <- bg[bg$gene %in% pc_symbols, , drop = FALSE]
cat("Background genes retained (protein-coding):", nrow(bg_pc), "\n")
cat("Background genes removed:", nrow(bg) - nrow(bg_pc), "\n")

# Write output
out_file <- file.path(script_dir, "background_genes_protein_coding.csv")
write.csv(bg_pc, out_file, row.names = FALSE)
cat("Output written to:", out_file, "\n")

sessionInfo()

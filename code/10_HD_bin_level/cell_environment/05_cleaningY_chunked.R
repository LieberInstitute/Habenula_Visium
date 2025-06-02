library(jaffelab)
library(here)
library(tidyverse)
library(scran)
library(spatialLIBD)
library(sessioninfo)
library(HDF5Array)
library(Matrix)

#   Get k from array task ID
k = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

#spe <- readRDS("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/ficture_harmony/spe_raw.rds")
spe <- readRDS("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/probe_fix/ficture_harmony/spe_raw.rds")

original_row_names <- rownames(spe)
original_col_names <- colnames(spe)
nonzero_rows <- rowSums(assays(spe)$counts) > 0
nonzero_cols <- colSums(assays(spe)$counts) > 0
#  Choose the nonzero bins
spe = spe[nonzero_rows, nonzero_cols]
#   Did a normalization before cleaningY()
#   Use library-size normalization (normalization by deconvolution is not
#   computationally feasible with data this large). Don't log scale, as for
#   FICTURE we want counts that statistically resemble real counts
message(Sys.time(), ' | Performing log normalization...')
spe = computeLibraryFactors(spe)
spe = logNormCounts(spe, transform = "none")
assays(spe)$counts) = assays(spe)$normcounts
assays(spe)$normcounts = NULL

mod <- with(colData(spe), model.matrix(~ sample_id))

# split the genes
n_genes <- nrow(spe)
gene_indices <- split(seq_len(n_genes), cut(seq_len(n_genes), 1000, labels = FALSE))

result_list <- vector("list", length = 1000)

for (i in ((1:100)+100*(k-1))) {
  cat("Cleaning chunk", i, "\n")
  
  idx <- gene_indices[[i]]
  counts_chunk <- assays(spe)$counts[idx, , drop = FALSE]
  zero_mask <- counts_chunk == 0
  # regression
  cleaned <- cleaningY(counts_chunk, mod, P = 1)
  cleaned[zero_mask] <- 0
  cleaned_sparse <- Matrix(cleaned, sparse = TRUE)
  result_list[[i]] <- cleaned_sparse

  rm(counts_chunk, cleaned, cleaned_sparse)
  gc()
}

final_result <- do.call(rbind, result_list)

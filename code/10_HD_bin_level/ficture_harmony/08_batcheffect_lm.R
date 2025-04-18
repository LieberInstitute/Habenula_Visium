#sample the genes 

library(jaffelab)
library(here)
library(tidyverse)
library(scran)
library(spatialLIBD)
library(sessioninfo)
library(HDF5Array)
library(Matrix)

spe <- readRDS("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/ficture_harmony/spe_raw.rds")
mod <- with(colData(spe), model.matrix(~ sample_id))

# split the genes
n_genes <- nrow(spe)
gene_indices <- split(seq_len(n_genes), cut(seq_len(n_genes), 1000, labels = FALSE))

# saving address
# dir.create("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/ficture_harmony/cleaned_chunks")

result_list <- vector("list", length = 1000)

for (i in seq_along(gene_indices)) {
  cat("Cleaning chunk", i, "\n")
  
  idx <- gene_indices[[i]]
  counts_chunk <- assays(spe)$counts[idx, , drop = FALSE]
  
  # regression
  cleaned <- cleaningY(counts_chunk, mod, P = 2)
  cleaned_sparse <- Matrix(cleaned, sparse = TRUE)
  result_list[[i]] <- cleaned_sparse

  rm(counts_chunk, cleaned, cleaned_sparse)
  gc()
}

final_result <- do.call(rbind, result_list)
assays(spe)$counts <- final_result

saveRDS(spe, file = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/ficture_harmony/spe/y_clean_spe.rds")

# #combine all
# cleaned_list <- lapply(list.files("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/ficture_harmony/cleaned_chunks", full.names = TRUE), readRDS)
# y_clean_p2 <- do.call(rbind, cleaned_list)
# dim(y_clean_p2)
# saveRDS(y_clean_p2, "y_clean_p2.rds")

# n_genes <- nrow(spe)
# gene_indices <- split(seq_len(n_genes), cut(seq_len(n_genes), 1000, labels = FALSE))
# spe_chunks <- lapply(gene_indices, function(idx) spe[idx, ])

# # Sample the genes
# n_genes <- nrow(spe)
# gene_indices <- split(seq_len(n_genes), cut(seq_len(n_genes), 1000, labels = FALSE))
# result_list <- vector("list", length = 1000)

# for (i in seq_along(gene_indices)) {
#   cat("Processing chunk", i, "\n")
  
#   idx <- gene_indices[[i]]
#   chunk <- spe[idx, ]  # 动态生成子集（不会创建整个chunks列表）
  
#   expr_mean <- rowMeans(assay(chunk, "counts"))
#   result_list[[i]] <- data.frame(
#     gene = rownames(chunk),
#     mean_count = expr_mean
#   )
  
#   rm(chunk)
#   gc()
# }

# final_result <- do.call(rbind, result_list)


# # Define a full model
# mod <- with(colData(spe), model.matrix(~ sample_id))
# y_clean_p2 <- cleaningY(assays(spe)$counts, mod, P=2)

# saveRDS(y_clean_p2, file = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/ficture_harmony/spe/y_clean_spe.rds")

# 
library(jaffelab)
library(here)
library(tidyverse)
library(scran)
library(spatialLIBD)
library(sessioninfo)
library(HDF5Array)
library(Matrix)

library(purrr)

k_vals <- 1:10

prefix <- "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/probe_fix/ficture_harmony/spe/count_cleaned_"

spe_list <- map(k_vals, ~ readRDS(paste0(prefix, .x, ".rds")))
spe_merged <- do.call(rbind, spe_list)

spe1 <- readRDS("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/probe_fix/ficture_harmony/spe_raw.rds")
original_row_names <- rownames(spe1)
original_col_names <- colnames(spe1)
nonzero_rows <- which(rowSums(assays(spe1)$counts) > 0)
nonzero_cols <- which(colSums(assays(spe1)$counts) > 0)

coo <- summary(spe_merged)
full_counts <- sparseMatrix(
  i = nonzero_rows[coo$i],
  j = nonzero_cols[coo$j],
  x = coo$x,
  dims = c(length(original_row_names), length(original_col_names)),
  dimnames = list(original_row_names, original_col_names)
)

assays(spe1)$counts <- full_counts

# range(full_counts)
# [1] 0.000000 8.566133

#saveRDS(spe, file = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/ficture_harmony/spe/y_clean_spe.rds")
saveRDS(spe1, file = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/probe_fix/ficture_harmony/spe/y_clean_spe.rds")

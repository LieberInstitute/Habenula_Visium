library(sessioninfo)
library(here)
library(SpatialExperiment)
library(scran)
library(scater)
library(tidyverse)

spe_in_path = here(
    'processed-data', '10_HD_bin_level', 'ficture_harmony', 'spe_raw.rds'
)
out_path = here(
    'processed-data', '10_HD_bin_level', 'ficture_harmony',
    'normalized_input.tsv.gz'
)

spe <- readRDS(spe_in_path)

#   Filter raw SPE: drop bins with 0 counts for all genes, and drop genes with
#   0 counts in every bin
message(Sys.time(), ' | Filtering bins and genes...')
spe <- spe[rowSums(assays(spe)$counts) > 0, colSums(assays(spe)$counts) > 0]

#   Use library-size normalization (normalization by deconvolution is not
#   computationally feasible with data this large). Don't log scale, as for
#   FICTURE we want counts that statistically resemble real counts
message(Sys.time(), ' | Performing log normalization...')
spe = computeLibraryFactors(spe)
spe = logNormCounts(spe, transform = "none")

assays(spe)$counts = NULL
gc()

#   Form a tibble of the nonzero elements of the normalized-counts matrix
#   (in a memory-efficient way)
message(Sys.time(), ' | Converting counts to FICTURE input...')
counts_mat = as(assays(spe)$logcounts, "TsparseMatrix")
counts_df = tibble(
        X = spatialCoords(spe)[counts_mat@j + 1, 1],
        Y = spatialCoords(spe)[counts_mat@j + 1, 2],
        gene = rownames(spe)[counts_mat@i + 1],
        Count = counts_mat@x,
        key = spe$key[counts_mat@j + 1]
    ) |>
    arrange(X)


write_tsv(counts_df, out_path)

session_info()

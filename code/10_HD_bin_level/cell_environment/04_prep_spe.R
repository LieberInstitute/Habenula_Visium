library(here)
library(scran)
library(tidyverse)
library(SpatialExperiment)
library(sessioninfo)

spe_in_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'ficture_harmony',
    'spe_raw.rds'
)
spe_out_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'cell_environment',
    'spe_filtered.rds'
)
bin_set_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'cell_environment',
    'extracellular_bins.csv.gz'
)

################################################################################
#   Take only extracellular bins and genes with nonzero expression
################################################################################

message(Sys.time(), ' | Loading...')
spe = readRDS(spe_in_path)
spe$key = paste(spe$sample_id, colnames(spe), sep = '_')

bin_set = read_csv(bin_set_path, show_col_types = FALSE) |>
    mutate(key = paste(sample_id, bin_id, sep = '_')) |>
    pull(key)

message(
    Sys.time(),
    sprintf(
        ' | Keeping extracellular bins (%.1f%%)',
        100 * mean(spe$key %in% bin_set)
    )
)
spe = spe[, spe$key %in% bin_set]

#   Filter raw SPE: drop bins with 0 counts for all genes, and drop genes with
#   0 counts in every bin
message(Sys.time(), ' | Dropping zero-expression bins and genes...')
spe = spe[rowSums(assays(spe)$counts) > 0, colSums(assays(spe)$counts) > 0]

################################################################################
#   Library-size normalize
################################################################################

#   Use library-size normalization (normalization by deconvolution is not
#   computationally feasible with data this large)
message(Sys.time(), ' | Performing log normalization...')
spe = computeLibraryFactors(spe)
spe = logNormCounts(spe)
assays(spe)$counts = NULL

message(Sys.time(), ' | Saving filtered SPE...')
saveRDS(spe, spe_out_path)

session_info()

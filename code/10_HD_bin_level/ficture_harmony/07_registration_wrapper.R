#   Deprecated script: begin spatial registration of cell-level FICTURE results.
#   We later determined registration could be directly performed with 2um
#   bin-level results, which is what we opted with going forward

library(here)
library(spatialLIBD)
library(HDF5Array)
library(sessioninfo)
library(tidyverse)

spe_dir = here('processed-data', '09_HD_cell_level', 'spe_norm_filtered')
cluster_path = here(
    'processed-data', '10_HD_bin_level', 'ficture_harmony',
    'ficture_merged.csv'
)
pseudo_path = here(
    'processed-data', '10_HD_bin_level', 'ficture_harmony',
    'pseudobulk_spe', 'library_normalized.rds'
)
model_path = here(
    'processed-data', '10_HD_bin_level', 'ficture_harmony',
    'modeling_results', 'library_normalized.rds'
)

dir.create(dirname(pseudo_path), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(model_path), showWarnings = FALSE)

#   Load and bring counts into memory to speed up computations. Despite the huge
#   size of the data, the memory footprint is manageable due to the extreme
#   sparsity of the data
spe = loadHDF5SummarizedExperiment(spe_dir)
assays(spe)$counts = as(assays(spe)$counts, "dgCMatrix")

#   Add in cluster assignments to 'spe', removing NA cells
cluster_df = read_csv(cluster_path, show_col_types = FALSE)
stopifnot(all(spe$key %in% cluster_df$key))
spe$FICTURE_cluster = cluster_df$FICTURE_cluster[match(spe$key, cluster_df$key)]
spe$FICTURE_cluster = factor(
    spe$FICTURE_cluster, levels = sort(unique(spe$FICTURE_cluster))
)
spe = spe[, !is.na(spe$FICTURE_cluster)]

#   Pseudobulk
model_results = registration_wrapper(
    spe,
    var_registration = 'FICTURE_cluster',
    var_sample_id = 'sample_id',
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name',
    pseudobulk_rds_file = pseudo_path
)

saveRDS(model_results, model_path)

session_info()

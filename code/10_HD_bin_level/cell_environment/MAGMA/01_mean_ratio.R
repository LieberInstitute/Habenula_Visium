#   Find mean-ratio markers for k = 10 FICTURE clusters. These will be input
#   gene sets for MAGMA

library(SpatialExperiment)
library(tidyverse)
library(DeconvoBuddies)
library(here)
library(sessioninfo)
library(spatialLIBD)
library(duckplyr)

k = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
spe_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'spe_filtered.rds'
)
cluster_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'ficture_plotting', 'extracellular.parquet'
)
out_dir = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'MAGMA', 'gene_sets'
)

mean_ratio_threshold = 1.05
max_num_genes = 200

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

spe = readRDS(spe_path)

ficture_df = read_parquet_duckdb(cluster_path, prudence = 'stingy') |>
    dplyr::rename(factor_K1 = paste0('k', k)) |>
    select(bin_key, factor_K1)

spe$ficture_cluster = tibble(
        bin_key = paste(colnames(spe), spe$sample_id, sep = "_"),
        idx = seq_along(colnames(spe))
    ) |>
    left_join(ficture_df, by = 'bin_key') |>
    #   This is critical, as duckplyr does not naturally preserve row order
    arrange(idx) |>
    pull(factor_K1)

rowData(spe)$gene_id = rownames(spe)
rowData(spe)$gene_name = rowData(spe)$symbol

spe = spe[, !is.na(spe$ficture_cluster)]

export_set = function(spe, cell_type_col, file_tag) {
    #   Pseudobulking accomplishes two things:
    #       1. Circumvents expression filtering in mean ratio which is likely
    #          too aggressive for this highly sparse data
    #       2. Applies basic expression filtering instead
    spe_pb = registration_pseudobulk(
        spe, var_registration = cell_type_col, var_sample_id = "sample_id"
    )
    
    marker_stats = get_mean_ratio(
            sce = spe_pb, assay_name = "logcounts",
            cellType_col = cell_type_col, gene_ensembl = "gene_id",
            gene_name = "gene_name"
        ) |>
        filter(MeanRatio > mean_ratio_threshold) |>
        dplyr::rename(set_id = cellType.target, gene_id = gene_ensembl) |>
        group_by(set_id) |>
        arrange(desc(MeanRatio)) |>
        slice_head(n = max_num_genes) |>
        ungroup() |>
        select(set_id, gene_id, gene_name, MeanRatio) |>
        arrange(set_id, desc(MeanRatio))

    message(sprintf("Marker counts for %s resolution:", file_tag))
    print(table(marker_stats$set_id))

    write_tsv(marker_stats, file.path(out_dir, sprintf("%s.tsv", file_tag)))
}

export_set(spe, "ficture_cluster", sprintf("k%d", k))

session_info()

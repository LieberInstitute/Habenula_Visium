# Find MeanRatio marker genes
library(SpatialExperiment)
library(tidyverse)
library(DeconvoBuddies)
library(here)
library(sessioninfo)
library(spatialLIBD)

spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'spe_norm_filtered_split.rds'
)
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'banksy', 'leiden_res1_8.csv'
)
cell_type_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy', 'cluster_annotation.csv'
)

out_dir = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA', 'gene_sets'
)

mean_ratio_threshold = 1.05
max_num_genes = 200

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

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
        dplyr::rename(set_id = cellType.target, gene_id = gene) |>
        group_by(set_id) |>
        arrange(desc(MeanRatio)) |>
        slice_head(n = max_num_genes) |>
        select(set_id, gene_id) |>
        arrange(set_id)

    message(sprintf("Marker counts for %s resolution:", file_tag))
    print(table(marker_stats$set_id))

    write_tsv(marker_stats, file.path(out_dir, sprintf("%s.tsv", file_tag)))
}

spe = readRDS(spe_path)

# ---- add banksy cluster -> cell type annotations into spe ----
cluster_df = read_csv(cluster_path, show_col_types = FALSE)
cell_type_df = read_csv(cell_type_path, show_col_types = FALSE)

cluster_df = cluster_df %>%
    left_join(cell_type_df, by = c("banksy" = "cluster"))

stopifnot(all(spe$key %in% cluster_df$key))

spe$cell_type_broad = cluster_df$broad_cell_type[match(spe$key, cluster_df$key)]
spe$cell_type_fine  = cluster_df$fine_cell_type [match(spe$key, cluster_df$key)]

spe$cell_type_broad = factor(spe$cell_type_broad, levels = sort(unique(spe$cell_type_broad)))
spe$cell_type_fine  = factor(spe$cell_type_fine,  levels = sort(unique(spe$cell_type_fine)))

# ---- export gene sets (fine / mid / broad) ----
export_set(spe, "cell_type_fine",  "fine")
export_set(spe, "cell_type_broad", "broad")

session_info()
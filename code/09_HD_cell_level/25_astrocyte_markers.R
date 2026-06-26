#   Are literature astrocyte markers expressed in the MHb?

library(here)
library(tidyverse)
library(spatialLIBD)
library(sessioninfo)

spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'no_secondary', 'astrocyte_markers'
)
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'banksy',
    'leiden_res1_8.csv'
)
ct_anno_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'cluster_annotation.csv'
)
sample_ids = c('Br9090_1', 'Br8433_1')
markers = c(
    "GFAP", "AQP4", "S100B", "SLC1A2", "SLC1A3", "APOE", "VIM", "NFIA", "NFIB"
)
cluster_colors = c(
    '6' = "#502419", '11' = '#19647E', '18' = '#F29559',
    'Other' = '#959595'
)

dir.create(plot_dir, showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

vis_clus_hd = function(
        spe, clustervar, sample_id, plot_path, colors = cluster_colors
    ) {
    #   Run twice to overcome a bug with different behavior on the first plot
    for (i in seq_len(2)) {
        p = vis_clus(
                spe, sampleid = sample_id, clustervar = clustervar,
                is_stitched = TRUE, point_size = 30, spatial = FALSE,
                colors = colors
            ) +
            guides(fill = guide_legend(override.aes = list(size = 8)))
    }
    png(plot_path, width = 1500, height = 1500)
    print(p)
    dev.off()
}

################################################################################
#   Main
################################################################################

spe = readRDS(spe_path)

stopifnot(all(markers %in% rowData(spe)$gene_name))
markers_ensembl = rowData(spe)$gene_id[match(markers, rowData(spe)$gene_name)]

anno_df = tibble(key = spe$key) |>
    left_join(read_csv(cluster_path, show_col_types = FALSE), by = 'key') |>
    dplyr::rename(cluster = banksy) |>
    left_join(read_csv(ct_anno_path, show_col_types = FALSE), by = 'cluster')
spe$banksy_cluster = ifelse(
    anno_df$cluster %in% c(6, 11, 18), as.character(anno_df$cluster), 'Other'
)
spe$cell_type = anno_df$fine_cell_type
stopifnot(!any(is.na(spe$cell_type)))

for (sample_id in sample_ids) {
    vis_clus_hd(
        spe, 'banksy_cluster', sample_id,
        file.path(plot_dir, sprintf('Astro_banksy_%s.png', sample_id))
    )

    p = vis_gene(
        spe, sampleid = sample_id, geneid = markers_ensembl,
        is_stitched = TRUE, point_size = 30, spatial = FALSE,
        cap_percentile = 0.99
    )
    png(
        file.path(plot_dir, sprintf('Astro_markers_%s.png', sample_id)),
        width = 1500, height = 1500
    )
    print(p)
    dev.off()
}

session_info()

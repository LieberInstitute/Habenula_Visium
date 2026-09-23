library(here)
library(tidyverse)
library(spatialLIBD)
library(sessioninfo)

sample_id = 'Br9090_1'
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'banksy',
    'leiden_res1_8.csv'
)
ct_anno_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'cluster_annotation.csv'
)
spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
plot_path = here(
    'plots', '09_HD_cell_level', 'no_secondary',
    sprintf('hb_only_banksy_%s.png', sample_id)
)
cell_type_colors = c(
    'MHb.1' = '#413C58',
    'MHb.2' = '#C98CA7',
    'LHb.1.3.4' = '#D5CB0A',
    'LHb.2.7' = '#C26D0D',
    'LHb.4' = '#9B1D20',
    'other' = '#C4C4C4'
)

################################################################################
#   Functions
################################################################################

vis_banksy = function(spe, clustervar, sample_id, plot_path) {
    #   Run twice to overcome a bug with different behavior on the first plot
    for (i in seq_len(2)) {
        p = vis_clus(
                spe, sampleid = sample_id, clustervar = clustervar,
                is_stitched = TRUE, point_size = 20, spatial = FALSE,
                colors = cell_type_colors
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

anno_df = read_csv(ct_anno_path, show_col_types = FALSE)

spe$cell_type = tibble(key = spe$key) |>
    left_join(read_csv(cluster_path, show_col_types = FALSE), by = 'key') |>
    mutate(
        cell_type = anno_df$fine_cell_type[
            match(as.character(banksy), anno_df$cluster)
        ]
    ) |>
    mutate(
        cell_type = ifelse(
                cell_type %in% names(cell_type_colors), cell_type, 'other'
            ) |>
            factor(levels = names(cell_type_colors))
    ) |>
    pull(cell_type)
stopifnot(!any(is.na(spe$cell_type)))

vis_banksy(
    spe = spe, clustervar = 'cell_type', sample_id = sample_id,
    plot_path = plot_path
)

session_info()

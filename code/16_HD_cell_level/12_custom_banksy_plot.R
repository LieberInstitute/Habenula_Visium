#   The Visium HD spatial registration data was used to disambiguate annotations
#   of some multiome clusters. This script plots some particular Banksy
#   resolutions and sets of clusters used in this disambiguation process. It
#   also plots some custom annotated plots for a talk for Kristen

library(here)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(sessioninfo)

sample_id_1 = "H1-MVPY9BW_A1_8433"
sample_id_2 = "H1-W369TJK_D1_9090"
spe_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'spe_norm_filtered.rds'
)
cluster_1_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'banksy', 'lambda0_2',
    'leiden_res1_2.csv'
)
cluster_2_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'banksy', 'lambda0_2',
    'leiden_res1_7.csv'
)
plot_1_path = here(
    'plots', '09_HD_cell_level', 'new_samples2', 'banksy', 'lambda0_2',
    'leiden_res1_2', sprintf('clusters_%s_custom.png', sample_id_1)
)
plot_2_path = here(
    'plots', '09_HD_cell_level', 'new_samples2', 'banksy', 'lambda0_2',
    'leiden_res1_7', sprintf('clusters_%s_hb_only.png', sample_id_2)
)
ct_anno_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'registration_banksy',
    'cluster_annotation.csv'
)
ambig_colors = c(
    "5" = "#35B42B",
    "13" = "#2F97FF",
    "14" = "#583E23",
    "18" = "#FFA239",
    "other" = "#DFE1DD"
)
anno_colors = c(
    MHb.1 = '#f38021',
    MHb.2 = '#3EA4FD',
    LHb.1.3.4 = '#004F2D',
    LHb.2.7 = '#00a900',
    Other = "#DFE1DD"
)

################################################################################
#   Functions
################################################################################

vis_banksy = function(spe, clustervar, colors, plot_path) {
    #   Run twice to overcome a bug with different behavior on the first plot
    for (i in seq_len(2)) {
        p = vis_clus(
                spe, clustervar = clustervar, is_stitched = TRUE,
                point_size = 20, spatial = FALSE, colors = colors
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

#   Mapping from clusters to cell types
anno_df = read_csv(ct_anno_path, show_col_types = FALSE)

#   Merge in Banksy clusters to SPE
temp = colnames(spe)
colData(spe) = colData(spe) |>
    as_tibble() |>
    left_join(
        read_csv(cluster_1_path, show_col_types = FALSE) |>
            dplyr::rename(banksy = banksy_lambda0_2),
        by = 'key'
    ) |>
    left_join(
        read_csv(cluster_2_path, show_col_types = FALSE) |>
            dplyr::rename(banksy2 = banksy_lambda0_2),
        by = 'key'
    ) |>
    mutate(
        banksy_ambig = factor(
            ifelse(
                as.character(banksy) %in% names(ambig_colors),
                as.character(banksy),
                "other"
            ),
            levels = names(ambig_colors)
        ),
        cell_type = anno_df$fine_cell_type[
            match(as.character(banksy2), anno_df$cluster)
        ],
        cell_type_subset = factor(
            ifelse(
                cell_type %in% names(anno_colors), cell_type, 'Other'
            ),
            levels = names(anno_colors)
        )
    ) |>
    DataFrame()
colnames(spe) = temp
stopifnot(!any(is.na(spe$banksy)))
stopifnot(!any(is.na(spe$banksy2)))

vis_banksy(
    spe[, spe$sample_id == sample_id_1], 'banksy_ambig', ambig_colors,
    plot_1_path
)

vis_banksy(
    spe[, spe$sample_id == sample_id_2], 'cell_type_subset', anno_colors,
    plot_2_path
)

session_info()

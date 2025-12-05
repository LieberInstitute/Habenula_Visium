#   The Visium HD spatial registration data was used to disambiguate annotations
#   of some multiome clusters. This script plots some particular Banksy
#   resolutions and sets of clusters used in this disambiguation process

library(here)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(sessioninfo)

sample_id = "H1-MVPY9BW_A1_8433"
spe_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'spe_norm_filtered.rds'
)
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'banksy', 'lambda0_2',
    'leiden_res1_2.csv'
)
plot_path = here(
    'plots', '09_HD_cell_level', 'new_samples2', 'banksy', 'lambda0_2',
    'leiden_res1_2', sprintf('clusters_%s_custom.png', sample_id)
)
cluster_colors = c(
    "5" = "#35B42B",
    "13" = "#2F97FF",
    "14" = "#583E23",
    "18" = "#FFA239",
    "other" = "#DFE1DD"
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

#   Load just the sample we're using to disambiguate
spe = readRDS(spe_path)
spe = spe[, spe$sample_id == sample_id]

#   Merge in Banksy clusters to SPE
temp = colnames(spe)
colData(spe) = colData(spe) |>
    as_tibble() |>
    left_join(
        read_csv(cluster_path, show_col_types = FALSE) |>
            dplyr::rename(banksy = banksy_lambda0_2),
        by = 'key'
    ) |>
    mutate(
        banksy = factor(
            ifelse(
                as.character(banksy) %in% names(cluster_colors),
                as.character(banksy),
                "other"
            ),
            levels = names(cluster_colors)
        )
    ) |>
    DataFrame()
colnames(spe) = temp
stopifnot(!any(is.na(spe$banksy)))

vis_banksy(spe, 'banksy', cluster_colors, plot_path)

session_info()

#   The Visium HD spatial registration data was used to disambiguate annotations
#   of some multiome clusters. This script plots some particular Banksy
#   resolutions and sets of clusters used in this disambiguation process

library(here)
library(tidyverse)
library(SpatialExperiment)
library(HDF5Array)
library(spatialLIBD)
library(sessioninfo)

sample_id = "H1-MVPY9BW_A1_8433"
spe_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples', 'spe_norm_filtered.rds'
)
cluster_1_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples', 'banksy', 'lambda0_2',
    'leiden_res1_4_subset.csv'
)
cluster_3_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples', 'banksy', 'lambda0_2',
    'leiden_res1_7_subset.csv'
)
plot_1_path = here(
    'plots', '09_HD_cell_level', 'new_samples', 'banksy', 'lambda0_2',
    'leiden_res1_4', sprintf('clusters_%s_subset_custom_1.png', sample_id)
)
plot_2_path = here(
    'plots', '09_HD_cell_level', 'new_samples', 'banksy', 'lambda0_2',
    'leiden_res1_4', sprintf('clusters_%s_subset_custom_2.png', sample_id)
)
plot_3_path = here(
    'plots', '09_HD_cell_level', 'new_samples', 'banksy', 'lambda0_2',
    'leiden_res1_7', sprintf('clusters_%s_subset_custom_1.png', sample_id)
)
cluster_colors_1 = c(
    "2" = "#35B42B",
    "3" = "#2F97FF",
    "9" = "#583E23",
    "11" = "#FFA239",
    "other" = "#DFE1DD"
)
cluster_colors_2 = c("17" = "#001DAF", "other" = "#DFE1DD")
cluster_colors_3 = c("5" = "#001DAF", "other" = "#DFE1DD")

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
spe = loadHDF5SummarizedExperiment(spe_dir)
spe = spe[, spe$sample_id == sample_id]

#   Merge in Banksy clusters to SPE
temp = colnames(spe)
colData(spe) = colData(spe) |>
    as_tibble() |>
    left_join(
        read_csv(cluster_1_path, show_col_types = FALSE) |>
            dplyr::rename(banksy_1_4 = banksy_lambda0_2),
        by = 'key'
    ) |>
    left_join(
        read_csv(cluster_3_path, show_col_types = FALSE) |>
            dplyr::rename(banksy_1_7 = banksy_lambda0_2),
        by = 'key'
    ) |>
    mutate(
        banksy_1 = factor(
            ifelse(
                as.character(banksy_1_4) %in% names(cluster_colors_1),
                as.character(banksy_1_4),
                "other"
            ),
            levels = names(cluster_colors_1)
        ),
        banksy_2 = factor(
            ifelse(
                as.character(banksy_1_4) %in% names(cluster_colors_2),
                as.character(banksy_1_4),
                "other"
            ),
            levels = names(cluster_colors_2)
        ),
        banksy_3 = factor(
            ifelse(
                as.character(banksy_1_7) %in% names(cluster_colors_3),
                as.character(banksy_1_7),
                "other"
            ),
            levels = names(cluster_colors_3)
        )
    ) |>
    DataFrame()
colnames(spe) = temp
stopifnot(!any(is.na(spe$banksy_1)))
stopifnot(!any(is.na(spe$banksy_2)))
stopifnot(!any(is.na(spe$banksy_3)))

#   Plot each custom plot
vis_banksy(spe, 'banksy_1', cluster_colors_1, plot_1_path)
vis_banksy(spe, 'banksy_2', cluster_colors_2, plot_2_path)
vis_banksy(spe, 'banksy_3', cluster_colors_3, plot_3_path)

session_info()

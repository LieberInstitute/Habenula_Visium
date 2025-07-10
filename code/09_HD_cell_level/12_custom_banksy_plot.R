#   The Visium HD spatial registration data was used to disambiguate annotations
#   of some multiome clusters. This script plots a particular Banksy resolution
#   and set of clusters used in this disambiguation process

library(here)
library(tidyverse)
library(SpatialExperiment)
library(HDF5Array)
library(spatialLIBD)
library(sessioninfo)

sample_id = "H1-MVPY9BW_A1_8433"
spe_dir = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'spe_norm_filtered'
)
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'banksy', 'lambda0_2',
    'leiden_res1_4_subset.csv'
)
plot_path = here(
    'plots', '09_HD_cell_level', 'probe_fix', 'banksy', 'lambda0_2',
    'leiden_res1_4', sprintf('clusters_%s_subset_custom.png', sample_id)
)
cluster_colors = c(
    "2" = "#246eb9",
    "3" = "#2BA822",
    "9" = "#583E23",
    "11" = "#ED8105",
    "other" = "#ABAFA9"
)

#   Load just the sample we're using to disambiguate
spe = loadHDF5SummarizedExperiment(spe_dir)
spe = spe[, spe$sample_id == sample_id]

#   Merge in Banksy clusters to SPE
spe$banksy = tibble(key = spe$key) |>
    left_join(read_csv(cluster_path, show_col_types = FALSE), by = 'key') |>
    mutate(
        banksy = ifelse(
            as.character(banksy_lambda0_2) %in% names(cluster_colors),
            as.character(banksy_lambda0_2),
            "other"
        )
    ) |>
    pull(banksy)
stopifnot(!any(is.na(spe$banksy)))

p = vis_clus(
        spe, clustervar = 'banksy', is_stitched = TRUE, point_size = 20,
        spatial = FALSE, colors = cluster_colors
    ) +
    guides(fill = guide_legend(override.aes = list(size = 8)))
png(plot_path, width = 1500, height = 1500)
print(p)
dev.off()

session_info()

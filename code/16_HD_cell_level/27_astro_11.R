#   Plot Astro cluster 11 (which located in MHb) spatially for a supplemental
#   figure

library(here)
library(tidyverse)
library(spatialLIBD)
library(sessioninfo)

cluster_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'banksy',
    'leiden_res1_8.csv'
)
spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'no_secondary', 'misc_paper_figs',
    'astro_11'
)
cell_map_path = here('raw-data', 'sample_info', 'hd_cell_type_map.csv')
sample_ids = c('Br3942_1', 'Br8433_2', 'Br8667_1', 'Br9090_2', 'Br9902_1')

dir.create(plot_dir, showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

vis_clus_hd = function(
        spe, clustervar, sample_id, plot_path, colors
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

spe$cell_type = tibble(key = spe$key) |>
    left_join(read_csv(cluster_path, show_col_types = FALSE), by = 'key') |>
    mutate(cell_type = ifelse(banksy == 11, '11 ~ Astrocyte', 'Other')) |>
    pull(cell_type)

astro_colors = c(
    '11 ~ Astrocyte' = read_csv(cell_map_path, show_col_types = FALSE) |>
        filter(new_cell_type == 'Astrocyte') |>
        pull(color),
    Other = '#9e9e9e'
)

for (this_sample_id in sample_ids) {
    vis_clus_hd(
        spe, clustervar = 'cell_type', sample_id = this_sample_id,
        plot_path = file.path(plot_dir, paste0(this_sample_id, '.png')),
        colors = astro_colors
    )
}

session_info()

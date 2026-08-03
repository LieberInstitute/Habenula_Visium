#   Plots for the supplemental figure HD_final_anno

library(here)
library(tidyverse)
library(spatialLIBD)
library(sessioninfo)

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
plot_dir = here(
    'plots', '09_HD_cell_level', 'no_secondary', 'misc_paper_figs',
    'final_anno_fig'
)
cell_map_path = here('raw-data', 'sample_info', 'hd_cell_type_map.csv')
sample_ids = c('Br8433_1', 'Br9090_1')

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

cell_df = tibble(key = spe$key) |>
    left_join(read_csv(cluster_path, show_col_types = FALSE), by = 'key') |>
    left_join(
        read_csv(ct_anno_path, show_col_types = FALSE),
        by = c('banksy' = 'cluster')
    ) |>
    left_join(
        read_csv(cell_map_path, show_col_types = FALSE),
        by = c('fine_cell_type' = 'old_cell_type')
    )
spe$cell_type = cell_df$new_cell_type
spe$banksy = cell_df$banksy
spe = spe[, !is.na(spe$cell_type)]

for (this_sample_id in sample_ids) {
    colors = c(
        '13' = '#BB2908', '14' = '#301BCB', 'Ependymal' = '#A5753C',
        'Other' = '#9e9e9e'
    )
    spe$temp = case_when(
            spe$cell_type == 'Ependymal' ~ 'Ependymal',
            spe$banksy %in% c(13, 14) ~ as.character(spe$banksy),
            TRUE ~ 'Other'
        ) |>
        factor(levels = names(colors))
    vis_clus_hd(
        spe, clustervar = 'temp', sample_id = this_sample_id,
        plot_path = file.path(plot_dir, paste0(this_sample_id, '_13_14.png')),
        colors = colors
    )

    colors = c(
        '9' = '#D006D7', 'Excit.Thal' = '#ED7F17', 'Other' = '#9e9e9e'
    )
    spe$temp = case_when(
            spe$cell_type == 'Excit.Thal' ~ 'Excit.Thal',
            spe$banksy == 9 ~ '9',
            TRUE ~ 'Other'
        ) |>
        factor(levels = names(colors))
    vis_clus_hd(
        spe, clustervar = 'temp', sample_id = this_sample_id,
        plot_path = file.path(plot_dir, paste0(this_sample_id, '_9.png')),
        colors = colors
    )
}

session_info()

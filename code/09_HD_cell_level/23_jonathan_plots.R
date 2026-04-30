#   Plot requests for Jonathan

library(here)
library(tidyverse)
library(spatialLIBD)
library(sessioninfo)
library(cowplot)

cluster_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'banksy',
    'leiden_res1_8.csv'
)
spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
anno_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'cluster_annotation.csv'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'no_secondary', 'jonathan_plots'
)
other_color = '#d4d4d4'
cluster_combos = list(
    MHb = c(
        MHb.1 = "#A86A9A",
        MHb.2 = "#BCA6B6",
        LHb = "#1f78b4",
        Other = other_color
    ),
    Excit_LHb = c(
        MHb = "#ad1d8c",
        LHb = "#1f78b4",
        Excit_LHb = "#A8B8BC",
        Other = other_color
    ),
    LHb.2.7 = c(
        MHb = "#ad1d8c",
        LHb = "#1f78b4",
        LHb.2.7 = "#6C9FA9",
        Other = other_color
    ),
    LHb.4 = c(
        MHb = "#ad1d8c",
        LHb = "#1f78b4",
        LHb.4 = "#00607A",
        Other = other_color
    ),
    Inhib_LHb_4.2 = c(
        MHb = "#ad1d8c",
        LHb = "#1f78b4",
        Inhib_LHb_4.2 = "#DC143C",
        Excit.Thal = "#4d55b7",
        Other = other_color
    ),
    Ependymal = c(
        Ependymal = "#f5a105",
        Subependymal = "#a5aa04",
        Other = other_color
    ),
    summary_plot = c(
        MHb = "#ad1d8c",
        LHb = "#1f78b4",
        Inhib_LHb_4.2 = "#DC143C",
        Excit.Thal = "#4d55b7",
        Other = other_color
    ),
    all_Hb = c(
        MHb = '#4024bb',
        Excit_LHb = '#D5CB0A',
        LHb.2.7 = '#C26D0D',
        LHb.4 = '#9B1D20',
        Inhib_LHb_4.2 = '#037d1f',
        Other = '#C4C4C4'
    )
)

dir.create(plot_dir, showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

#   vis_clus with much less white space
vis_clus_improved = function(spe, sampleid, clustervar, colors) {
    small_spe = spe[,spe$sample_id == sampleid]

    p = spatialCoords(small_spe) |>
        as_tibble() |>
        mutate(
            x = pxl_col_in_fullres,
            y = max(pxl_row_in_fullres) - pxl_row_in_fullres,
            cluster = small_spe[[clustervar]]
        ) |>
        ggplot(aes(x = x, y = y, color = cluster)) +
            geom_point(size = 0.3) +
            scale_color_manual(values = colors) +
            coord_fixed() +
            labs(color = 'Cell Type', title = sampleid) +
            guides(color = guide_legend(override.aes = list(size = 5))) +
            theme_void(base_size = 15)
    
    return(p)
}

plot_combo = function(spe, combo_name) {
    stopifnot(!is.null(spe[[combo_name]]))

    spe[[combo_name]] = factor(
        spe[[combo_name]], levels = names(cluster_combos[[combo_name]])
    )

    all_samples = unique(spe$sample_id)[grepl('_1$', unique(spe$sample_id))]

    p_list = list()
    for (this_sample_id in all_samples) {
        #   Run twice to overcome a bug with different behavior on the first
        #   plot
        for (i in seq_len(2)) {
            p_list[[this_sample_id]] = vis_clus_improved(
                    spe, sampleid = this_sample_id, clustervar = combo_name,
                    colors = cluster_combos[[combo_name]]
                )
        }
    }
    p = plot_grid(plotlist = p_list, nrow = 1)
    png(
        file.path(plot_dir, sprintf('%s.png', combo_name)),
        width = 2000, height = 400
    )
    print(p)
    dev.off()
}

################################################################################
#   Main
################################################################################

spe = readRDS(spe_path)

anno_df = read_csv(anno_path, show_col_types = FALSE)
spe$banksy = tibble(key = spe$key) |>
    left_join(read_csv(cluster_path, show_col_types = FALSE), by = 'key') |>
    pull(banksy)
stopifnot(!any(is.na(spe$banksy)))

#-------------------------------------------------------------------------------
#   Define combinations of clusters to view based on Jonathan's request
#-------------------------------------------------------------------------------

spe$MHb = case_when(
    spe$banksy %in% c(7, 19, 23, 25) ~ "LHb",
    spe$banksy == 8 ~ "MHb.1",
    spe$banksy == 15 ~ "MHb.2",
    TRUE ~ "Other"
)
spe$Excit_LHb = case_when(
    spe$banksy %in% c(7, 19, 25) ~ "LHb",
    spe$banksy %in% c(8, 15) ~ "MHb",
    spe$banksy == 23 ~ "Excit_LHb",
    TRUE ~ "Other"
)
spe$LHb.2.7 = case_when(
    spe$banksy %in% c(7, 23, 25) ~ "LHb",
    spe$banksy %in% c(8, 15) ~ "MHb",
    spe$banksy == 19 ~ "LHb.2.7",
    TRUE ~ "Other"
)
spe$LHb.4 = case_when(
    spe$banksy %in% c(19, 23) ~ "LHb",
    spe$banksy %in% c(8, 15) ~ "MHb",
    spe$banksy %in% c(7, 25) ~ "LHb.4",
    TRUE ~ "Other"
)
spe$Inhib_LHb_4.2 = case_when(
    spe$banksy %in% c(7, 19, 23, 25) ~ "LHb",
    spe$banksy %in% c(8, 15) ~ "MHb",
    spe$banksy == 20 ~ "Excit.Thal",
    spe$banksy %in% c(5, 9) ~ "Inhib_LHb_4.2",
    TRUE ~ "Other"
)
spe$Ependymal = case_when(
    spe$banksy %in% c(12, 24) ~ "Ependymal",
    spe$banksy %in% c(13, 14) ~ "Subependymal",
    TRUE ~ "Other"
)
spe$summary_plot = case_when(
    spe$banksy %in% c(7, 19, 23, 25) ~ "LHb",
    spe$banksy %in% c(8, 15) ~ "MHb",
    spe$banksy == 20 ~ "Excit.Thal",
    spe$banksy %in% c(5, 9) ~ "Inhib_LHb_4.2",
    TRUE ~ "Other"
)
spe$all_Hb = case_when(
    spe$banksy == 23 ~ "Excit_LHb",
    spe$banksy == 19 ~ "LHb.2.7",
    spe$banksy %in% c(7, 25) ~ "LHb.4",
    spe$banksy %in% c(8, 15) ~ "MHb",
    spe$banksy %in% c(5, 9) ~ "Inhib_LHb_4.2",
    TRUE ~ "Other"
)

#-------------------------------------------------------------------------------
#   Plot combinations
#-------------------------------------------------------------------------------

for (combo_name in names(cluster_combos)) {
    plot_combo(spe, combo_name)
}

session_info()

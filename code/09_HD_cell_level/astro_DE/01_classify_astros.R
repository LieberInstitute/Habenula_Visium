library(tidyverse)
library(here)
library(jsonlite)
library(RANN)
library(sessioninfo)

spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'banksy',
    'leiden_res1_8.csv'
)
ct_anno_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'cluster_annotation.csv'
)
out_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'astro_DE',
    'astro_labels.csv.gz'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'no_secondary', 'astro_DE',
    'astro_classification'
)
example_samples = c("Br9090_1", "Br8433_1")
sample_info_path = here('raw-data', 'sample_info', 'hd_basic_info_split.csv')
ASTRO_DIST_THRESHOLD_UM = 100
ASTRO_KNN = 10
label_colors = c(
    "Astro: medial" = "#2DF2FD",
    "Astro: lateral" = "#EAD637",
    "Astro: neither" = "#000000",
    "MHb" = "#1E3888",
    "LHb" = "#BE0625",
    "Other" = "#979797"
)

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_path), showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

vis_clus_custom = function(spe, clustervar, sample_id, plot_path) {
    #   Run twice to overcome a bug with different behavior on the first plot
    for (i in seq_len(2)) {
        p = vis_clus(
                spe, sampleid = sample_id, clustervar = clustervar,
                is_stitched = TRUE, point_size = 20, spatial = FALSE,
                colors = label_colors
            ) +
            guides(fill = guide_legend(override.aes = list(size = 8)))
    }
    png(plot_path, width = 1500, height = 1500)
    print(p)
    dev.off()
}

################################################################################
#   Join in cell-type annotation
################################################################################

spe = readRDS(spe_path)

spe$cell_type = tibble(key = spe$key) |>
    left_join(read_csv(cluster_path, show_col_types = FALSE), by = 'key') |>
    dplyr::rename(cluster = banksy) |>
    left_join(read_csv(ct_anno_path, show_col_types = FALSE), by = 'cluster') |>
    pull(fine_cell_type)
stopifnot(!any(is.na(spe$cell_type)))

################################################################################
#   Convert spatial coordinates into units of microns
################################################################################

sample_info_df = read_csv(sample_info_path, show_col_types = FALSE) |>
    select(tissue_id, spaceranger_dir)

#   Build per-sample microns_per_pixel lookup from spaceranger JSON
mpp_lookup = sample_info_df |>
    mutate(
        json_path = here(
            spaceranger_dir,
            "outs/binned_outputs/square_002um/spatial/scalefactors_json.json"
        ),
        microns_per_pixel = map_dbl(
            json_path, \(p) fromJSON(p)$microns_per_pixel
        )
    ) |>
    select(tissue_id, microns_per_pixel)

#   Convert spatial coordinates from pixels to microns, per sample
mpp_per_spot = mpp_lookup$microns_per_pixel[
    match(spe$sample_id, mpp_lookup$tissue_id)
]

coord_df = cbind(
        as.data.frame(spatialCoords(spe) * mpp_per_spot),
        sample_id = spe$sample_id,
        cell_type = spe$cell_type,
        key = spe$key
    ) |>
    as_tibble() |>
    filter(grepl('^Astrocyte$|^MHb|LHb', cell_type)) |>
    mutate(
        cell_type = case_when(
            cell_type == 'Astrocyte' ~ 'Astrocyte',
            grepl('^MHb', cell_type) ~ 'MHb',
            grepl('LHb', cell_type) ~ 'LHb',
            TRUE ~ NA_character_
        )
    )
stopifnot(!any(is.na(coord_df$cell_type)))

################################################################################
#   Classify astrocytes as medial, lateral, or neither
################################################################################

#   For each sample, compute mean distance from each astrocyte to its nearest
#   10 MHb and 10 LHb neurons using a KD-tree, then threshold at 100 µm.
#   Astrocytes within threshold of both nuclei are assigned to whichever is
#   closer on average.

COORDS = c("pxl_col_in_fullres", "pxl_row_in_fullres")

astro_labels = purrr::map(
        unique(coord_df$sample_id),
        function(s) {
            df_s = dplyr::filter(coord_df, sample_id == s)

            astro_m = as.matrix(df_s[df_s$cell_type == "Astrocyte", COORDS])
            mhb_m   = as.matrix(df_s[df_s$cell_type == "MHb",       COORDS])
            lhb_m   = as.matrix(df_s[df_s$cell_type == "LHb",       COORDS])

            tibble(
                key           = df_s$key[df_s$cell_type == "Astrocyte"],
                mean_dist_mhb = rowMeans(
                    nn2(mhb_m, astro_m, k = ASTRO_KNN)$nn.dists
                ),
                mean_dist_lhb = rowMeans(
                    nn2(lhb_m, astro_m, k = ASTRO_KNN)$nn.dists
                )
            )
        }
    ) |>
    list_rbind()

astro_df = coord_df |>
    filter(cell_type == "Astrocyte") |>
    left_join(astro_labels, by = "key") |>
    mutate(
        astro_label = case_when(
            mean_dist_mhb < ASTRO_DIST_THRESHOLD_UM &
                mean_dist_lhb < ASTRO_DIST_THRESHOLD_UM ~
                    ifelse(mean_dist_mhb < mean_dist_lhb, "medial", "lateral"),
            mean_dist_mhb < ASTRO_DIST_THRESHOLD_UM ~ "medial",
            mean_dist_lhb < ASTRO_DIST_THRESHOLD_UM ~ "lateral",
            TRUE ~ "neither"
        )
    ) |>
    select(key, astro_label)

################################################################################
#   Plots visually validating that our classification is reasonable
################################################################################

#   Choice of distance threshold
p = astro_labels |>
    pivot_longer(
        c(mean_dist_mhb, mean_dist_lhb),
        names_to = "nucleus", values_to = "mean_dist"
    ) |>
    mutate(nucleus = if_else(nucleus == "mean_dist_mhb", "MHb", "LHb")) |>
    ggplot(aes(mean_dist, color = nucleus, fill = nucleus)) +
        geom_density(alpha = 0.2) +
        geom_vline(xintercept = ASTRO_DIST_THRESHOLD_UM, linetype = "dashed") +
        scale_x_continuous(
            #   Cap at 99th percentile of LHb (the tighter distribution) so
            #   both peaks are readable without the long MHb tail dominating
            limits = c(0, quantile(astro_labels$mean_dist_lhb, 0.99))
        ) +
        labs(
            x = sprintf("Mean distance to nearest %d neurons (µm)", ASTRO_KNN),
            y = "Density", color = NULL, fill = NULL
        ) +
        theme_bw(base_size = 15)
pdf(file.path(plot_dir, "astro_distance_density.pdf"))
print(p)
dev.off()

spe$astro_label = tibble(key = spe$key) |>
    left_join(astro_df, by = "key") |>
    mutate(
        astro_label = case_when(
            is.na(astro_label) & grepl('^MHb', spe$cell_type) ~ "MHb",
            is.na(astro_label) & grepl('LHb', spe$cell_type) ~ "LHb",
            is.na(astro_label) ~ "Other",
            TRUE ~ paste('Astro:', astro_label)
        )
    ) |>
    pull(astro_label)

for (this_sample in example_samples) {
    vis_clus_custom(
        spe = spe, clustervar = 'astro_label', sample_id = this_sample,
        plot_path = file.path(
            plot_dir, sprintf('astro_labels_%s.png', this_sample)
        )
    )
}
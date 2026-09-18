library(here)
library(tidyverse)
library(crawdad)
library(spatialLIBD)
library(scales)
library(sessioninfo)

sample_info_path = here('raw-data', 'sample_info', 'hd_basic_info_split.csv')
result_paths = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'crawdad',
    'output', '%s_results.csv'
)
in_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'crawdad',
    'input_cells.csv.gz'
)
spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
plot_dir = here('plots', '09_HD_cell_level', 'no_secondary', 'crawdad')
cell_type_colors = c(
    Astrocyte = "#2F97FF", 'LHb.2.7' = "#FFA239", Other = "#DFE1DD"
)
min_num_signif = 5
cell_type_levels = c(
    'MHb.1', 'MHb.2', 'Excit_LHb', 'LHb.2.7', 'LHb.4', 'LHb.4/Inhib_LHb_4.2',
    'Inhib_LHb_4.2', 'Astrocyte', 'Endo', 'Endo/microglia', 'Oligo', 'OPC',
    'Ependymal', 'Subependymal'
)

dir.create(
    file.path(plot_dir, 'spatial_plots'), recursive = TRUE, showWarnings = FALSE
)

################################################################################
#   Functions
################################################################################

custom_dotplot = function(result_df, z_sig, filename) {
    p = ggplot(
            result_df, aes(x = reference, y = neighbor, color = Z, size = scale)
        ) +
        geom_point() +
        scale_color_gradientn(
            colors = c('blue', '#CECECE', '#CECECE', 'red'),
            values = rescale(
                c(min(result_df$Z), -1 * z_sig, z_sig, max(result_df$Z))
            )
        ) +
        scale_radius(
            trans = 'reverse',
            breaks = seq(
                min(result_df$scale), max(result_df$scale), length.out = 3
            ),
            range = c(2, 8)
        ) +
        coord_fixed() +
        theme_bw(base_size = 20) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    pdf(file.path(plot_dir, filename), width = 9)
    print(p)
    dev.off()
}

################################################################################
#   CRAWDAD-specific plots
################################################################################

sample_info = read_csv(sample_info_path, show_col_types = FALSE)
sample_ids = sample_info$tissue_id

result_list = list()
for (sample_id in sample_ids) {
    result_list[[sample_id]] = sprintf(result_paths, sample_id) |>
        read_csv(show_col_types = FALSE) |>
        mutate(sample_id = sample_id)
}

#   Get the Z-score significance threshold (same in all samples)
z_sig = do.call(rbind, result_list) |>
    filter(sample_id == sample_ids[1]) |>
    correctZBonferroni()

result_df = do.call(rbind, result_list) |>
    filter(reference != 'Excit.Thal', neighbor != 'Excit.Thal') |>
    mutate(
        reference = factor(reference, levels = cell_type_levels),
        neighbor = factor(neighbor, levels = cell_type_levels)
    ) |>
    #   First average Z-scores across permutations
    group_by(sample_id, neighbor, scale, reference) |>
    summarize(Z = mean(Z)) |>
    ungroup() |>
    #   Then filter to the smallest spatial scale with significant Z-scores
    filter(abs(Z) >= z_sig) |>
    group_by(sample_id, neighbor, reference) |>
    filter(scale == min(scale)) |>
    #   Retain pairs where all samples all signs of Z scores agree across
    #   samples, and significance is achieved in some sufficient number of
    #   samples 
    group_by(reference, neighbor) |>
    filter(all(Z > 0) | all(Z < 0)) |>
    filter(n() >= min_num_signif) |>
    #   Take the mean Z-score and scale across samples
    group_by(neighbor, reference) |>
    summarize(scale = mean(scale), Z = mean(Z)) |>
    ungroup() |>
    #   Cap Z-score at twice the magnitude of the significance threshold
    mutate(Z = sign(Z) * pmin(abs(Z), z_sig * 2))

custom_dotplot(result_df, z_sig, 'dot_plot.pdf')

#   Plot Z-scores vs scale for a particularly interesting cell-type pair
p = do.call(rbind, result_list) |>
    #   Average Z-scores across permutations
    group_by(sample_id, neighbor, scale, reference) |>
    summarize(Z = mean(Z)) |>
    ungroup() |>
    #   Improve plot appearance
    mutate(
        facet_anno = sprintf("Ref: %s\nNeighbor: %s", reference, neighbor)
    ) |>
    #   Focus on a particular pair (and its reverse)
    filter(
        ((neighbor == 'Astrocyte') & (reference == 'LHb.2.7')) |
        ((neighbor == 'LHb.2.7') & (reference == 'Astrocyte'))
    ) |>
    ggplot(aes(x = scale, y = Z, color = sample_id, group = sample_id)) +
        geom_line() +
        geom_point() +
        geom_hline(yintercept = z_sig, linetype = 'dashed') +
        geom_hline(yintercept = -1 * z_sig, linetype = 'dashed') +
        facet_wrap(~ facet_anno, nrow = 1) +
        theme_bw(base_size = 20) +
        labs(x = 'Scale (Microns)', color = 'Sample ID')

pdf(
    file.path(plot_dir, 'z_scores_Astrocyte_LHb_2_7.pdf'),
    width = 10, height = 5
)
print(p)
dev.off()

cell_df = read_csv(in_path, show_col_types = FALSE) |>
    mutate(
        cell_type = ifelse(
            cell_type %in% c('Astrocyte', 'LHb.2.7'), cell_type, 'Other'
        )
    ) |>
    select(key, cell_type)

spe = readRDS(spe_path)
spe = spe[, spe$key %in% cell_df$key]

spe$cell_type = tibble(key = spe$key) |>
    left_join(cell_df, by = 'key') |>
    pull(cell_type)

#   Plot the cell-type pair spatially in each sample
dir.create(
    file.path(plot_dir, 'spatial_plots', 'Astrocyte_LHb_2_7'),
    showWarnings = FALSE
)
for (sample_id in sample_ids) {
    #   Run twice to overcome a bug with different behavior on the first
    #   plot
    for (i in seq_len(2)) {
        p = vis_clus(
                spe, sampleid = sample_id, clustervar = 'cell_type',
                is_stitched = TRUE, point_size = 20, spatial = FALSE,
                colors = cell_type_colors
            ) +
            guides(fill = guide_legend(override.aes = list(size = 8)))
    }
    png(
        file.path(
            plot_dir, 'spatial_plots', 'Astrocyte_LHb_2_7',
            sprintf('%s.png', sample_id)
        ),
        width = 1500, height = 1500
    )
    print(p)
    dev.off()
}

session_info()

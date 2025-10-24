library(here)
library(tidyverse)
library(crawdad)
library(spatialLIBD)
library(scales)
library(sessioninfo)

sample_info_path = here('raw-data', 'sample_info', 'hd_basic_info.csv')
result_paths = here(
    'processed-data', '09_HD_cell_level', 'new_samples', 'crawdad', 'region',
    'output', '%s_%s_results.csv'
)
in_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples', 'crawdad', 'region',
    'input_cells.csv.gz'
)
spe_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples', 'spe_norm_filtered.rds'
)
plot_dir = here('plots', '09_HD_cell_level', 'new_samples', 'crawdad', 'region')
regions = c('habenula', 'thalamus')
cell_type_colors = c(
    Astrocyte = "#2F97FF", placeholder = "#FFA239", Other = "#DFE1DD"
)

dir.create(file.path(plot_dir, 'spatial_plots'), showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

custom_dotplot = function(result_df, z_sig, cell_types, filename) {
    p = ggplot(
            result_df, aes(x = reference, y = neighbor, color = Z, size = scale)
        ) +
        geom_point() +
        scale_x_discrete(limits = cell_types) +
        scale_y_discrete(limits = cell_types) +
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
            range = c(2, 15)
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

sample_info = read_csv(sample_info_path)
sample_ids = sample_info$sample_id[1:3]

result_list = list()
for (sample_id in sample_ids) {
    for (region in regions) {
        result_list[[paste0(sample_id, region)]] = sprintf(
                result_paths, sample_id, region
            ) |>
            read_csv(show_col_types = FALSE) |>
            mutate(sample_id = sample_id, region = region)
    }
}

#   Get the Z-score significance threshold (same in all samples/regions)
z_sig = do.call(rbind, result_list) |>
    filter(region == regions[1], sample_id == sample_ids[1]) |>
    correctZBonferroni()

result_df = do.call(rbind, result_list) |>
    #   First average Z-scores across permutations
    group_by(region, sample_id, neighbor, scale, reference) |>
    summarize(Z = mean(Z)) |>
    ungroup() |>
    #   Then filter to the smallest spatial scale with significant Z-scores
    filter(abs(Z) >= z_sig) |>
    group_by(region, sample_id, neighbor, reference) |>
    filter(scale == min(scale)) |>
    #   Retain pairs where all samples are significant, and sign of Z scores
    #   agree across samples
    group_by(region, reference, neighbor) |>
    filter(all(Z > 0) | all(Z < 0)) |>
    filter(n() == length(sample_ids)) |>
    #   Take the mean Z-score and scale across samples
    group_by(region, neighbor, reference) |>
    summarize(scale = mean(scale), Z = mean(Z)) |>
    ungroup() |>
    #   Cap Z-score at twice the magnitude of the significance threshold
    mutate(Z = sign(Z) * pmin(abs(Z), z_sig * 2))

#   Custom dot plot for each region
for (region_name in regions) {
    #   Grab all unique cell types originally present in the data
    cell_types = do.call(rbind, result_list) |>
        filter(region == region_name) |>
        pull(reference) |>
        unique() |>
        sort()

    result_df |>
        filter(region == region_name) |>
        custom_dotplot(
            z_sig, cell_types, sprintf('dot_plot_%s.pdf', region_name)
        )
}

for (hb_subtype in c('MHb.2', 'LHb.7')) {
    #   Plot Z-scores vs scale for a particularly interesting cell-type pair
    p = do.call(rbind, result_list) |>
        #   Average Z-scores across permutations
        group_by(region, sample_id, neighbor, scale, reference) |>
        summarize(Z = mean(Z)) |>
        ungroup() |>
        #   Improve plot appearance
        mutate(
            sample_id = paste0('Br', str_extract(sample_id, '[0-9]{4}$')),
            facet_anno = sprintf("Ref: %s\nNeighbor: %s", reference, neighbor)
        ) |>
        #   Focus on a particular pair (and its reverse)
        filter(
            region == 'habenula',
            ((neighbor == hb_subtype) & (reference == 'Astrocyte')) |
            ((neighbor == 'Astrocyte') & (reference == hb_subtype))
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
        file.path(plot_dir, sprintf('z_scores_%s_astro.pdf', hb_subtype)),
        width = 10, height = 5
    )
    print(p)
    dev.off()

    cell_df = read_csv(in_path, show_col_types = FALSE) |>
        filter(region_anno == 'habenula') |>
        mutate(
            cell_type = ifelse(
                cell_type %in% c(hb_subtype, 'Astrocyte'), cell_type, 'Other'
            )
        ) |>
        select(key, cell_type)

    spe = readRDS(spe_path)
    spe = spe[, spe$key %in% cell_df$key]

    spe$cell_type = tibble(key = spe$key) |>
        left_join(cell_df, by = 'key') |>
        pull(cell_type)

    #   Plot the cell-type pair spatially (only habenula) in each sample
    names(cell_type_colors) = c('Astrocyte', hb_subtype, 'Other')
    dir.create(
        file.path(plot_dir, 'spatial_plots', sprintf('%s_astro', hb_subtype)),
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
                plot_dir, 'spatial_plots', sprintf('%s_astro', hb_subtype),
                sprintf('%s.png', sample_id)
            ),
            width = 1500, height = 1500
        )
        print(p)
        dev.off()
    }
}

session_info()

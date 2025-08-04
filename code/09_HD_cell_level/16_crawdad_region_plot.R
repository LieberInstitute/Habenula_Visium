library(here)
library(tidyverse)
library(crawdad)
library(spatialLIBD)
library(HDF5Array)
library(scales)
library(sessioninfo)

sample_id_path = here('raw-data', 'sample_info', 'hd_sample_list.txt')
result_paths = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'crawdad', 'region',
    'output', '%s_%s_results.csv'
)
plot_dir = here('plots', '09_HD_cell_level', 'probe_fix', 'crawdad', 'region')
regions = c('habenula', 'thalamus')

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

sample_ids = readLines(sample_id_path)[1:3]

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

#   Grab all unique cell types for later
cell_types = do.call(rbind, result_list) |>
    pull(reference) |>
    unique() |>
    sort()

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
    result_df |>
        filter(region == region_name) |>
        custom_dotplot(
            z_sig, cell_types, sprintf('dot_plot_%s.pdf', region_name)
        )
}

session_info()

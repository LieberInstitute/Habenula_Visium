library(here)
library(tidyverse)
library(crawdad)
library(spatialLIBD)
library(scales)
library(sessioninfo)

k = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

sample_info_path = here('raw-data', 'sample_info', 'hd_basic_info_split.csv')
result_paths = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'crawdad', 'output', sprintf('%s_k%d_results.csv', '%s', k)
)
in_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'crawdad', sprintf('input_cells_k%d.csv.gz', k)
)
spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
plot_dir = here(
    'plots', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'crawdad'
)
cell_map_path = here('raw-data', 'sample_info', 'hd_cell_type_map.csv')
min_num_signif = 4

dir.create(plot_dir, showWarnings = FALSE)

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
            range = c(1, 5)
        ) +
        coord_fixed() +
        scale_x_discrete(drop = FALSE) +
        scale_y_discrete(drop = FALSE) +
        theme_bw(base_size = 16) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    pdf(file.path(plot_dir, filename), width = 9)
    print(p)
    dev.off()
}

################################################################################
#   Main
################################################################################

sample_info = read_csv(sample_info_path, show_col_types = FALSE)
sample_ids = sample_info |>
    filter(donor != 'Br9902') |>
    pull(tissue_id)

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

cell_map_df = rbind(
    read_csv(cell_map_path, show_col_types = FALSE) |>
        select(old_cell_type, new_cell_type),
    tibble(
        old_cell_type = paste0('Factor_', 0:(k - 1)),
        new_cell_type = paste0('Factor_', 0:(k - 1))
    )
)

result_df = do.call(rbind, result_list) |>
    filter(reference != 'Excit.Thal', neighbor != 'Excit.Thal') |>
    #   Use latest cell-type names and order properly
    mutate(
        reference = factor(
            cell_map_df$new_cell_type[
                match(reference, cell_map_df$old_cell_type)
            ],
            levels = cell_map_df$new_cell_type
        ),
        neighbor = factor(
            cell_map_df$new_cell_type[
                match(neighbor, cell_map_df$old_cell_type)
            ],
            levels = cell_map_df$new_cell_type
        )
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

custom_dotplot(result_df, z_sig, sprintf('dot_plot_k%d.pdf', k))

#   Do a version with just Banksy-derived cell types
if (k == 3) {
    result_df |>
        filter(!grepl('^Factor', reference), !grepl('^Factor', neighbor)) |>
        mutate(
            reference = droplevels(reference),
            neighbor = droplevels(neighbor)
        ) |>
        custom_dotplot(z_sig, 'dot_plot_cell_types.pdf')
}

session_info()

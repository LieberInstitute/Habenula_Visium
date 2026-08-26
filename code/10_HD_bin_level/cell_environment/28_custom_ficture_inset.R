library(tidyverse)
library(here)
library(duckplyr)
library(sessioninfo)
library(Polychrome)

k = 17
set.seed(6751)

extra_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'ficture_plotting', 'extracellular.parquet'
)
extra_bin_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'extracellular_bins.csv.gz'
)
plot_dir = here(
    'plots', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'ficture_plotting', sprintf('k_%d', k), 'custom_inset'
)
cluster_levels = paste0('Factor_', seq(0, k - 1))
factor_colors = setNames(
    Polychrome::glasbey.colors(k + 1)[2:(k + 1)], cluster_levels
)

num_cores = as.integer(Sys.getenv('SLURM_CPUS_PER_TASK'))
if (is.na(num_cores)) {
    num_cores = parallel::detectCores(logical = FALSE)
}
duckplyr::db_exec(sprintf('SET threads = %d', num_cores))
fallback_config(info = FALSE)

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

ficture_df = read_parquet_duckdb(extra_path, prudence = 'lavish') |>
    mutate(sample_id = str_extract(bin_key, '_(H1-.*)$', group = 1)) |>
    select(x, y, sample_id, bin_key, !!paste0('k', k)) |>
    dplyr::rename(ficture_cluster = !!paste0('k', k)) |>
    mutate(ficture_cluster = paste0('Factor_', ficture_cluster))

bin_cell_df = read_csv_duckdb(extra_bin_path, prudence = 'lavish') |>
    mutate(
        sample_id = str_extract(cell_key, '_(H1-.*)$', group = 1),
        bin_key = paste(bin_id, sample_id, sep = '_')
    ) |>
    select(bin_key, cell_key)

# Sample one mapped bin, then only collect nearby bins from the same sample.
# This keeps the large bin-to-cell join restricted to a small spatial window.
anchor_bin = bin_cell_df |>
    slice_sample(n = 1) |>
    inner_join(ficture_df, by = 'bin_key') |>
    collect()

sample_id = anchor_bin$sample_id
center_x = anchor_bin$x
center_y = anchor_bin$y
window_width = 500
max_cells = 35
min_bins = 100

sample_bin_cell_df = bin_cell_df |>
    filter(str_detect(cell_key, paste0('_', sample_id, '$'))) |>
    select(bin_key, cell_key)

window_df = ficture_df |>
    filter(
        sample_id == !!sample_id,
        between(x, !!center_x - window_width / 2, !!center_x + window_width / 2),
        between(y, !!center_y - window_width / 2, !!center_y + window_width / 2)
    ) |>
    inner_join(sample_bin_cell_df, by = 'bin_key') |>
    collect()

if (nrow(window_df) < min_bins) {
    stop('Random inset did not contain enough mapped bins. Try another seed or increase window_width.')
}

inset_widths = seq(window_width, 100, by = -25)
inset_summaries = tibble(inset_width = inset_widths) |>
    rowwise() |>
    mutate(
        n_bins = sum(
            between(window_df$x, center_x - inset_width / 2, center_x + inset_width / 2) &
                between(window_df$y, center_y - inset_width / 2, center_y + inset_width / 2)
        ),
        n_cells = n_distinct(window_df$cell_key[
            between(window_df$x, center_x - inset_width / 2, center_x + inset_width / 2) &
                between(window_df$y, center_y - inset_width / 2, center_y + inset_width / 2)
        ])
    ) |>
    ungroup()

inset_width = inset_summaries |>
    filter(n_cells <= max_cells, n_bins >= min_bins) |>
    arrange(desc(inset_width)) |>
    pull(inset_width) |>
    first()

if (is.na(inset_width)) {
    inset_width = inset_summaries |>
        filter(n_bins >= min_bins) |>
        arrange(n_cells, desc(n_bins)) |>
        pull(inset_width) |>
        first()
}

inset_df = window_df |>
    filter(
        between(x, center_x - inset_width / 2, center_x + inset_width / 2),
        between(y, center_y - inset_width / 2, center_y + inset_width / 2)
    )

cell_levels = sort(unique(inset_df$cell_key))
cell_colors = setNames(
    Polychrome::glasbey.colors(length(cell_levels) + 1)[2:(length(cell_levels) + 1)],
    cell_levels
)

save_inset_plot = function(plot_df, color_var, colors, legend_title, plot_path) {
    p = ggplot(plot_df, aes(x = x, y = 0 - y, color = {{ color_var }})) +
        geom_point(size = 2, shape = 15) +
        scale_color_manual(values = colors) +
        coord_fixed() +
        theme_void(base_size = 15) +
        labs(color = legend_title) +
        guides(color = 'none')
    ggsave(plot_path, p, dpi = 300, bg = 'white')

    return(invisible(NULL))
}

save_inset_plot(
    inset_df,
    ficture_cluster,
    factor_colors,
    'FICTURE factor',
    file.path(plot_dir, sprintf('%s_ficture_cluster_inset.pdf', sample_id))
)

save_inset_plot(
    inset_df,
    cell_key,
    cell_colors,
    'Cell',
    file.path(plot_dir, sprintf('%s_cell_inset.pdf', sample_id))
)

session_info()

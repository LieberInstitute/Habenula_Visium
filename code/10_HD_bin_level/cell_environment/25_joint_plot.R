library(tidyverse)
library(here)
library(duckplyr)
library(sessioninfo)
library(Polychrome)

k = c(4, 10, 20)[as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))]

extra_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'ficture_plotting', 'extracellular.parquet'
)
cellular_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'ficture_plotting', 'cellular.parquet'
)
plot_dir = here(
    'plots', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'ficture_plotting', sprintf('k_%d', k)
)
cluster_levels = paste0('Factor_', seq(0, k - 1))
factor_colors = setNames(
    Polychrome::palette.36[seq_along(cluster_levels)],
    cluster_levels
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

extra_df = read_parquet_duckdb(extra_path) |>
    mutate(sample_id = str_extract(bin_key, '_(H1-.*)$', group = 1)) |>
    select(x, y, sample_id, !!paste0('k', k)) |>
    dplyr::rename(ficture_cluster = !!paste0('k', k)) |>
    mutate(ficture_cluster = paste0('Factor_', ficture_cluster)) |>
    collect()

for (sample_id in unique(extra_df$sample_id)) {
    p = extra_df |>
        filter(sample_id == !!sample_id) |>
        ggplot(aes(x = x, y = 0 - y, color = ficture_cluster)) +
            geom_point(size = 0.2, shape = 15) +
            scale_color_manual(values = factor_colors) +
            coord_fixed() +
            theme_void(base_size = 15) +
            labs(color = 'FICTURE factor') +
            guides(color = guide_legend(override.aes = list(size = 10)))
    ggsave(
        file.path(plot_dir, sprintf('%s.png', sample_id)),
        p, width = 8, height = 7, dpi = 500
    )
}

session_info()

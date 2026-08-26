library(tidyverse)
library(here)
library(duckplyr)
library(sessioninfo)
library(Polychrome)

ficture_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
    'k8_cluster_coords.parquet'
)
plot_dir = here(
    'plots', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
    'ficture_plotting'
)
k = 8
cluster_levels = paste0('Factor_', seq(0, k - 1))
factor_colors = setNames(
    Polychrome::glasbey.colors(k+1)[2:(k+1)], cluster_levels
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(plot_dir, showWarnings = FALSE)

save_ficture_plot = function(extra_df, plot_path, colors = factor_colors) {
    p = ggplot(extra_df, aes(x = x, y = 0 - y, color = ficture_cluster)) +
        geom_point(size = 0.01, shape = 15) +
        scale_color_manual(values = colors) +
        coord_fixed() +
        theme_void(base_size = 15) +
        labs(color = 'FICTURE factor') +
        guides(color = guide_legend(override.aes = list(size = 10)))
    ggsave(plot_path, p, width = 24, height = 20, dpi = 150, bg = 'white')
  
    return(invisible(NULL))
}

ficture_df = read_parquet_duckdb(ficture_path, prudence = 'lavish') |>
    mutate(ficture_cluster = paste0('Factor_', ficture_cluster)) |>
    collect()

for (sample_id in unique(ficture_df$sample_id)) {
    ficture_df |>
        filter(sample_id == !!sample_id) |>
        save_ficture_plot(file.path(plot_dir, sprintf('%s.png', sample_id)))
}

session_info()

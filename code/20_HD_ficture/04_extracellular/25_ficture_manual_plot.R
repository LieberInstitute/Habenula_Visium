library(tidyverse)
library(here)
library(duckplyr)
library(sessioninfo)
library(Polychrome)

k = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

extra_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'ficture_plotting', 'extracellular.parquet'
)
plot_dir = here(
    'plots', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'ficture_plotting', sprintf('k_%d', k)
)
cluster_levels = paste0('Factor_', seq(0, k - 1))
factor_colors = setNames(
    Polychrome::glasbey.colors(k+1)[2:(k+1)], cluster_levels
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

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

extra_df = read_parquet_duckdb(extra_path, prudence = 'lavish') |>
    mutate(sample_id = str_extract(bin_key, '_(H1-.*)$', group = 1)) |>
    select(x, y, sample_id, !!paste0('k', k)) |>
    dplyr::rename(ficture_cluster = !!paste0('k', k)) |>
    mutate(ficture_cluster = paste0('Factor_', ficture_cluster)) |>
    collect()

for (sample_id in unique(extra_df$sample_id)) {
    extra_df |>
        filter(sample_id == !!sample_id) |>
        save_ficture_plot(file.path(plot_dir, sprintf('%s.png', sample_id)))
  
    if (k == 17) {
        extra_df |>
            filter(sample_id == !!sample_id) |>
            mutate(
                ficture_cluster = ifelse(
                    grepl('^Factor_[0-4]$', ficture_cluster), ficture_cluster, 'Other'
                )
            ) |>
            save_ficture_plot(
                file.path(plot_dir, sprintf('%s_5factor.png', sample_id)),
                colors = c(factor_colors[1:5], Other = 'grey80')
            )
        
        dir.create(file.path(plot_dir, 'individual'), showWarnings = FALSE)
        for (this_cluster in paste0('Factor_', seq(0, 4))) {
            extra_df |>
                filter(sample_id == !!sample_id) |>
                mutate(
                    ficture_cluster = ifelse(
                        ficture_cluster == !!this_cluster, this_cluster, 'Other'
                    )
                ) |>
                save_ficture_plot(
                    file.path(
                        plot_dir, 'individual',
                        sprintf('%s_%s.png', sample_id, this_cluster)
                    ),
                    colors = c(
                        setNames(factor_colors[1], this_cluster),
                        Other = 'grey80'
                    )
                )
        }
    }
}

session_info()

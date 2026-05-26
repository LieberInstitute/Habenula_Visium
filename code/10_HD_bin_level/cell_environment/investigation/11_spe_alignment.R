#   As a likely final check for bugs in the code, I'll check if FICTURE clusters
#   align correctly with the extracellular SPE by plotting them (crucially,
#   using spatial coordinates from the SPE, which has not been checked yet)

library(here)
library(spatialLIBD)
library(sessioninfo)
library(tidyverse)
library(duckplyr)
library(Polychrome)

k = 10
this_sample_id = 'H1-W369TJK_D1_9090'
spe_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'spe_filtered.rds'
)
cluster_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'ficture_plotting', 'extracellular.parquet'
)
cluster_raw_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'ficture_outputs', 'cleaningy', 'k_%d', 'analysis', 'nF%d.d_12',
    'cleaningy_joined_input.tsv.gz'
) |> sprintf(k, k)
plot_dir = here(
    'plots', '10_HD_bin_level', 'no_secondary', 'cell_environment', 
    'investigation'
)
cluster_levels = paste0('Factor_', seq(0, k - 1))
factor_colors = setNames(
    Polychrome::palette36.colors(k + 2)[3:(k + 2)], cluster_levels
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

################################################################################
#   Functions
################################################################################

plot_ficture = function(spe, plot_path) {
    p = tibble(
            x = spatialCoords(spe)[, 'pxl_col_in_fullres'],
            y = spatialCoords(spe)[, 'pxl_row_in_fullres'],
            ficture_cluster = paste0('Factor_', spe$ficture_cluster),
            sample_id = spe$sample_id
        ) |>
        filter(sample_id == this_sample_id) |>
        ggplot(aes(x = x, y = 0 - y, color = ficture_cluster)) +
            geom_point(size = 0.01, shape = 15) +
            scale_color_manual(values = factor_colors) +
            coord_fixed() +
            theme_void(base_size = 15) +
            labs(color = 'FICTURE factor') +
            guides(color = guide_legend(override.aes = list(size = 10)))
    ggsave(
        plot_path, p, width = 24, height = 20, dpi = 150, bg = 'white'
    )
}

################################################################################
#   Main
################################################################################

spe = readRDS(spe_path)

ficture_df = read_parquet_duckdb(cluster_path, prudence = 'stingy') |>
    dplyr::rename(factor_K1 = paste0('k', k)) |>
    select(bin_key, factor_K1)

spe$ficture_cluster = tibble(
        bin_key = paste(colnames(spe), spe$sample_id, sep = "_")
    ) |>
    left_join(ficture_df, by = 'bin_key') |>
    pull(factor_K1)

message(
    sprintf(
        "Dropping %.1f%% of bins without a defined FICTURE cluster",
        mean(is.na(spe$ficture_cluster)) * 100
    )
)
spe = spe[, !is.na(spe$ficture_cluster)]

plot_ficture(spe, file.path(plot_dir, 'alignment_my_parquet.png'))

#   Wow, I confirmed the bug! Is something already problematic from the outputs
#   directly from spatula? Above I used a more thoroughly processed version

ficture_df = read_csv_duckdb(cluster_raw_path, prudence = 'stingy') |>
    distinct(barcode, sample_id, factor_K1) |>
    filter(factor_K1 != 'NA')

spe$ficture_cluster = tibble(
        sample_id = spe$sample_id, barcode = colnames(spe)
    ) |>
    left_join(ficture_df, by = c('sample_id', 'barcode')) |>
    pull(factor_K1)

plot_ficture(spe, file.path(plot_dir, 'alignment_spatula.png'))

#   Extremely weird... this is still wrong, but these are different results
#   For my sanity I'll check the joining without duckplyr given the results
#   so far

session_info()

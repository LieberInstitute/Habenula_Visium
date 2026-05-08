#   The tiny, random-looking correlations from 02_cor_heatmap are suspicious;
#   one possible explanation that would explain all weird results seen so far
#   is improper matching of 2um bins with cells. Use plotting to visually
#   assess this mapping

library(tidyverse)
library(here)
library(duckplyr)
library(sessioninfo)
library(Polychrome)

ficture_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'ficture_plotting', 'extracellular.parquet'
)
bin_cell_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'extracellular_bins.csv.gz'
)
banksy_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'banksy',
    'leiden_res1_8.csv'
)
ct_anno_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'cluster_annotation.csv'
)
plot_dir = here(
    'plots', '10_HD_bin_level', 'no_secondary', 'cell_environment', 
    'investigation'
)
cluster_levels = c(paste0('Factor_', seq(0, 9)), 'outside_MHb')
factor_colors = setNames(
    c(Polychrome::palette36.colors(12)[3:12], "gray"), cluster_levels
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

anno_df = read_csv_duckdb(ct_anno_path, prudence = 'stingy') |>
    dplyr::rename(banksy = cluster, cell_type = fine_cell_type)

cluster_df = read_csv_duckdb(banksy_path, prudence = 'stingy') |>
    dplyr::rename(cell_key = key) |>
    left_join(anno_df, by = 'banksy') |>
    select(cell_key, cell_type)

ficture_df = read_parquet_duckdb(ficture_path, prudence = 'stingy') |>
    dplyr::rename(ficture_cluster = k10) |>
    filter(!is.na(ficture_cluster)) |>
    select(bin_key, x, y, ficture_cluster)

full_df = read_csv_duckdb(bin_cell_path, prudence = 'lavish') |>
    mutate(
        sample_id = str_extract(cell_key, '_(H1-.*)$', group = 1),
        bin_key = paste(bin_id, sample_id, sep = '_')
    ) |>
    select(cell_key, bin_key) |>
    inner_join(ficture_df, by = 'bin_key') |>
    inner_join(cluster_df, by = 'cell_key') |>
    collect()

p = full_df |>
    mutate(
        ficture_cluster = ifelse(
            grepl('^MHb', cell_type), paste0('Factor_', ficture_cluster),
            'outside_MHb'
        )
    ) |>
    filter(grepl('9090$', cell_key)) |>
    ggplot(aes(x = x, y = 0 - y, color = ficture_cluster)) +
        geom_point(size = 0.01, shape = 15) +
        scale_color_manual(values = factor_colors) +
        coord_fixed() +
        theme_void(base_size = 15) +
        labs(color = 'FICTURE factor') +
        guides(color = guide_legend(override.aes = list(size = 10)))
ggsave(
    file.path(plot_dir, 'MHb_Br9090.png'),
    p, width = 24, height = 20, dpi = 150, bg = 'white'
)

session_info()

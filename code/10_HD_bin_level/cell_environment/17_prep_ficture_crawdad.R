#   We plan to run a version of CRAWDAD incorporating clusters both from Banksy
#   and FICTURE. This script generates the inputs for CRAWDAD with both types
#   of clusters

library(tidyverse)
library(here)
library(duckplyr)
library(sessioninfo)
library(SpatialExperiment)
library(rjson)

k = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

extra_bin_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'extracellular_bins.csv.gz'
)
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'banksy',
    'leiden_res1_8.csv'
)
sample_info_path = here('raw-data', 'sample_info', 'hd_basic_info_split.csv')
anno_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'cluster_annotation.csv'
)
spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
ficture_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'ficture_plotting', 'extracellular.parquet'
)
crawdad_out_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'crawdad', sprintf('input_cells_k%d.csv.gz', k)
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(dirname(crawdad_out_path), showWarnings = FALSE)

ficture_df = read_parquet_duckdb(ficture_path, prudence = 'stingy') |>
    dplyr::rename(factor_K1 = paste0('k', k))

#   Mapping from clusters to cell types
anno_df = read_csv(anno_path, show_col_types = FALSE)

input_df = read_csv_duckdb(extra_bin_path, prudence = 'lavish') |>
    mutate(
        sample_id = str_extract(cell_key, '_(H1-.*)$', group = 1),
        bin_key = paste(bin_id, sample_id, sep = '_')
    ) |>
    inner_join(ficture_df, by = 'bin_key') |>
    #   For each cell, take the most common extracellular FICTURE cluster,
    #   randomly breaking any ties
    summarize(
        factor_K1 = {
            counts <- table(factor_K1)
            candidates <- names(counts[counts == max(counts)])
            sample(candidates, 1)
        },
        .by = cell_key
    ) |>
    #   Grab Banksy clusters for each cell
    inner_join(
        read_csv_duckdb(cluster_path, prudence = 'stingy') |>
            dplyr::rename(cell_key = key),
        by = 'cell_key'
    ) |>
    filter(banksy != 16) |>
    mutate(
        cell_type = anno_df$fine_cell_type[
            match(as.character(banksy), anno_df$cluster)
        ]
    ) |>
    collect()

sample_info = read_csv(sample_info_path, show_col_types = FALSE)

#   Ultimately, we'll be converting spatial coordinates to units of microns,
#   which is more interpretable than pixels
micron_per_px = c()
for (tissue_id in sample_info$tissue_id) {
    spaceranger_dir = sample_info$spaceranger_dir[
        sample_info$tissue_id == tissue_id
    ]
    scalefactors_path = here(
        spaceranger_dir, 'outs', 'binned_outputs', 'square_002um', 'spatial',
        'scalefactors_json.json'
    )
    micron_per_px = c(
        micron_per_px,
        fromJSON(file = scalefactors_path)[['microns_per_pixel']]
    )
}
scale_df = tibble(
    sample_id = sample_info$tissue_id,
    micron_per_px = micron_per_px
)

message('Number of bins per donor:')
table(input_df$sample_id)

#   Gather spatial coordinates, FICTURE, and Banksy clusters into a single CSV
#   for input to CRAWDAD
input_df |>
    dplyr::rename(ficture_cluster = factor_K1) |>
    mutate(
        x = x * scale_df$micron_per_px[match(sample_id, sample_info$tissue_id)],
        y = y * scale_df$micron_per_px[match(sample_id, sample_info$tissue_id)]
    ) |>
    select(x, y, sample_id, cell_type, ficture_cluster) |>
    write_csv(crawdad_out_path)

session_info()

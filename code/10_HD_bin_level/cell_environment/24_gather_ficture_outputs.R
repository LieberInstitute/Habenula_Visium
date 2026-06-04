#   Default FICTURE outputs have one row per gene per bin, which is a lot of
#   extra data, and slow to work with. Join all FICTURE outputs into a single
#   parquet with bin ID, coordinates, and cluster assignments for each k.
#   Next, export a parquet with one row per cellular bin, with bin ID,
#   coordinates, and cell type.
#
#   Overall, this script serves two functions:
#       1. Make FICTURE outputs easier to work with
#       2. Prepare the necessary data for jointly plotting FICTURE and Banksy
#          results

library(tidyverse)
library(here)
library(duckplyr)
library(SpatialExperiment)
library(sessioninfo)

k_values = 3:20

ficture_paths = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'ficture_outputs', 'cleaningy', 'k_%d', 'analysis', 'nF%d.d_12',
    'cleaningy_joined_input.tsv.gz'
)
cellular_bins_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'cellular_bins.csv.gz'
)
spe_extra_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'spe_filtered.rds'
)
spe_cell_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
banksy_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'banksy',
    'leiden_res1_8.csv'
)
ct_anno_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'cluster_annotation.csv'
)
extra_out_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'ficture_plotting', 'extracellular.parquet'
)
cellular_out_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'ficture_plotting', 'cellular.parquet'
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(dirname(extra_out_path), showWarnings = FALSE)

################################################################################
#   Extracellular bins
################################################################################

#   Columns 'bin_key', 'x', 'y', 'k[number]', ..., 'k[number]' for each k in
#   k_values

ficture_list = vector("list", length(k_values))
names(ficture_list) = paste0("k", k_values)

for (this_k in k_values) {
    ficture_list[[paste0("k", this_k)]] = sprintf(
            ficture_paths, this_k, this_k
        ) |>
        read_csv_duckdb(prudence = 'lavish') |>
        distinct(barcode, sample_id, factor_K1) |>
        filter(factor_K1 != 'NA') |>
        mutate(
            bin_key = paste(barcode, sample_id, sep = '_'),
            factor_K1 = as.integer(factor_K1)
        ) |>
        select(bin_key, !!paste0("k", this_k) := factor_K1) |>
        collect()
}

ficture_joined = purrr::reduce(ficture_list, full_join, by = "bin_key")

spe_extra = readRDS(spe_extra_path)

extra_df = tibble(
        bin_key = spe_extra$key,
        x = unname(spatialCoords(spe_extra)[, 1]),
        y = unname(spatialCoords(spe_extra)[, 2])
    ) |>
    #   Accidentally switched sample and bin ID when building the SPE
    mutate(
        bin_key = str_replace(bin_key, '^(H1.*_[0-9]{4})_(s_.*)$', '\\2_\\1')
    ) |>
    inner_join(ficture_joined, by = "bin_key") |>
    compute_parquet(extra_out_path)

################################################################################
#   Cellular bins
################################################################################

#   Columns 'bin_key', 'x', 'y', 'cell_type'

spe_cell = readRDS(spe_cell_path)

anno_df = read_csv(ct_anno_path, show_col_types = FALSE)

cell_df = tibble(
        key = spe_cell$key,
        x = unname(spatialCoords(spe_cell)[, 1]),
        y = unname(spatialCoords(spe_cell)[, 2])
    ) |>
    left_join(read_csv(banksy_path, show_col_types = FALSE), by = 'key') |>
    mutate(
        cell_type = anno_df$fine_cell_type[
            match(as.character(banksy), anno_df$cluster)
        ]
    ) |>
    dplyr::rename(cell_key = key) |>
    filter(cell_type != 'Drop') |>
    select(cell_key, x, y, cell_type)

bin_df = read_csv(cellular_bins_path, show_col_types = FALSE) |>
    dplyr::rename(bin_key = bin_id) |>
    inner_join(cell_df, by = "cell_key") |>
    select(bin_key, x, y, cell_type) |>
    compute_parquet(cellular_out_path)

session_info()

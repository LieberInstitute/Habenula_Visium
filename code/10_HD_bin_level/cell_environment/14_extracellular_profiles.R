library(tidyverse)
library(here)
library(spatialLIBD)
library(duckdb)
library(sessioninfo)

spe_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'spe_norm_filtered.rds'
)
extra_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'cell_environment',
    'extracellular_bins.csv.gz'
)
ficture_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
    'bin_level_clusters_batch.csv.gz'
)
banksy_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'banksy', 'lambda0_2',
    'leiden_res1_7.csv'
)
ct_anno_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'registration_banksy',
    'cluster_annotation.csv'
)
ficture_colnames = c('sample_id', 'barcode', 'FICTURE_k23')

spe = readRDS(spe_path)

#   Merge in annotated Banksy results
anno_df = read_csv(ct_anno_path, show_col_types = FALSE)
spe$cell_type = tibble(key = spe$key) |>
    left_join(read_csv(banksy_path, show_col_types = FALSE), by = 'key') |>
    mutate(
        cell_type = anno_df$fine_cell_type[
            match(as.character(banksy_lambda0_2), anno_df$cluster)
        ]
    ) |>
    pull(cell_type)

con = dbConnect(duckdb())

col_data = tibble(
    cell_id = spe$key, cell_type = spe$cell_type, bins_per_cell = spe$bin_count
)
duckdb_register(con, "col_data", col_data)

bin_df = fread(extra_path) |>
    as_tibble() |>
    dplyr::rename(barcode = bin_id) |>
    mutate(cell_id = paste(cell_id, sample_id, sep = '_')) |>
    left_join(col_data, by = 'cell_id') |>
    left_join(
        fread(ficture_path, select = ficture_colnames) |>
            as_tibble(),
        by = c('sample_id', 'barcode')
    )

sql_query = sprintf(
    "
    SELECT 
        extra.sample_id,
        extra.bin_id AS barcode,
        extra.cell_id || '_' || extra.sample_id AS cell_id,
        col_data.cell_type,
        col_data.bins_per_cell,
        ficture.FICTURE_k23
    FROM read_csv_auto('%s') AS extra
    LEFT JOIN col_data 
        ON extra.cell_id || '_' || extra.sample_id = col_data.cell_id
    LEFT JOIN read_csv_auto('%s') AS ficture
        ON extra.sample_id = ficture.sample_id AND extra.bin_id = ficture.barcode
    ",
    extra_path, ficture_path
)

bin_df = dbGetQuery(con, sql_query) |>
    as_tibble()

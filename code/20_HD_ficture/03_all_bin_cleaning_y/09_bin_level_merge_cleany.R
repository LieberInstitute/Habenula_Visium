#   Take all FICTURE results (across all k) and merge with the 2um bin-level
#   SPE, exporting a CSV of top clusters for each k (for easier handling than
#   the massive raw FICTURE results). Use duckdb for speed and memory

library(here)
library(tidyverse)
library(data.table)
library(SpatialExperiment)
library(sessioninfo)
library(duckdb)

ficture_input_paths = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
    'ficture_outputs', 'cleany', 'k_%d', 'analysis', 'nF%d.d_12',
    'transcripts_ficture_joined.tsv.gz'
)
spe_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony', 'spe',
    'y_clean_spe.rds'
)
out_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
    'bin_level_clusters_batch.csv.gz'
)
k_values = c(seq(3, 40), 70, 100)

ficture_colnames = c('sample_id', 'barcode', 'factor_K1')

con = dbConnect(duckdb())

#   Read in just the sample IDs and barcodes from the SPE
message(Sys.time(), ' | Grabbing and sorting colData...')
spe = readRDS(spe_path)
col_data = colData(spe) |>
    as_tibble() |>
    select(sample_id) |>
    mutate(barcode = colnames(spe))
rm(spe)
gc()
duckdb_register(con, "col_data", col_data)

#   Sort col_data to match with results from SQL queries
sorted_col_data = con |>
    dbGetQuery(
        "
        SELECT sample_id, barcode
        FROM col_data
        ORDER BY sample_id, barcode
        "
    ) |>
    as_tibble()

for (k in k_values) {
    message(Sys.time(), sprintf(' | Reading in k = %d results...', k))
    sql_query = sprintf(
        "
        SELECT ficture.factor_K1 AS cluster
        FROM col_data
        LEFT JOIN (
            SELECT DISTINCT sample_id, barcode, factor_K1
            FROM read_csv_auto('%s', delim = '\t')
            WHERE factor_K1 != 'NA'
        ) ficture
        ON col_data.sample_id = ficture.sample_id AND col_data.barcode = ficture.barcode
        ORDER BY col_data.sample_id, col_data.barcode
        ",
        sprintf(ficture_input_paths, k, k)
    )
    sorted_col_data[[paste0('FICTURE_k', k)]] = dbGetQuery(con, sql_query)$cluster
}

message(Sys.time(), ' | Writing to disk')
write_csv(sorted_col_data, out_path)

session_info()

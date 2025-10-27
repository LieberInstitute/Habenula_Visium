#   Take all FICTURE results (across all k) and merge with the 2um bin-level
#   SPE, exporting a CSV of top clusters for each k (for easier handling than
#   the massive raw FICTURE results)

library(here)
library(tidyverse)
library(data.table)
library(SpatialExperiment)
library(sessioninfo)

ficture_input_paths = here(
    'processed-data', '10_HD_bin_level', 'new_samples', 'ficture_harmony',
    'ficture_outputs', 'normalized', 'k_%d', 'analysis', 'nF%d.d_12',
    'normalized_joined_input.tsv.gz'
)
spe_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples', 'ficture_harmony',
    'spe_raw.rds'
)
out_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples', 'ficture_harmony',
    'bin_level_clusters.csv.gz'
)
k_values = c(seq(3, 40), 70, 100)

ficture_colnames = c('sample_id', 'barcode', 'factor_K1')

#   Read in just the sample IDs and barcodes from the SPE
spe = readRDS(spe_path)
col_data = colData(spe) |>
    as_tibble() |>
    select(sample_id) |>
    mutate(barcode = colnames(spe))
rm(spe)
gc()

for (k in k_values) {
    message(Sys.time(), sprintf(' | Reading in k = %d results...', k))
    ficture_df = sprintf(ficture_input_paths, k, k) |>
        #   Read in quickly but format as a tibble
        fread(select = ficture_colnames, sep = '\t') |>
        as_tibble() |>
        #   Set data types and drop empty rows
        mutate(sample_id = factor(sample_id), factor_K1 = factor(factor_K1)) |>
        filter(!is.na(factor_K1)) |>
        #   Name cluster column with value of k
        dplyr::rename(!!quo_name(paste0('FICTURE_k', k)) := factor_K1)
    
    col_data = left_join(
        col_data, ficture_df, by = c('sample_id', 'barcode'), multiple = 'any'
    )
    rm(ficture_df); gc()
}

message(Sys.time(), ' | Writing to disk')
write_csv(col_data, out_path)

session_info()

library(here)
library(tidyverse)
library(data.table)
library(SpatialExperiment)
library(sessioninfo)

k = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

ficture_path = here(
        'processed-data', '10_HD_bin_level', 'probe_fix', 'cell_environment',
        'ficture_outputs', 'normalized', 'k_%d', 'analysis', 'nF%d.d_12',
        'normalized_joined_input.tsv.gz'
    ) |>
    sprintf(k, k)
extra_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'cell_environment',
    'extracellular_bins.csv.gz'
)
out_path = here(
        'processed-data', '10_HD_bin_level', 'probe_fix', 'cell_environment',
        'cell_profiles', 'k_%d.csv.gz'
    ) |>
    sprintf(k)
ficture_colnames = c(
    'sample_id', 'barcode', 'factor_K1', 'factor_K2', 'factor_K3', 'factor_P1',
    'factor_P2', 'factor_P3'
)

dir.create(dirname(out_path), showWarnings = FALSE)

message(Sys.time(), ' | Reading in FICTURE clusters...')
ficture_df = fread(ficture_path, select = ficture_colnames) |>
    as_tibble() |>
    filter(!is.na(factor_K1))

message(Sys.time(), ' | Reading in extracellular bins and joining...')
extra_df = fread(extra_path) |>
    as_tibble() |>
    dplyr::rename(barcode = bin_id) |>
    #   This is much faster than first taking unique combintations of sample_id
    #   and barcode before joining
    inner_join(ficture_df, by = c('sample_id', 'barcode'), multiple = 'any')

#   Form score columns for each cluster based on adding up posterior
#   probabilities for the top 3 factors
for (this_k in seq_len(k) - 1) {
    extra_df[[paste0('score_', this_k)]] = (
        (extra_df$factor_K1 == this_k) * extra_df$factor_P1 +
        (extra_df$factor_K2 == this_k) * extra_df$factor_P2 +
        (extra_df$factor_K3 == this_k) * extra_df$factor_P3
    )
}

extra_df = extra_df |>
    #   For each cell and sample, sum up scores
    group_by(sample_id, cell_id) |>
    summarize(across(matches('^score_'), sum)) |>
    ungroup() |>
    #   Normalize scores so they add to 1 across all clusters
    mutate(temp_sum = rowSums(across(matches('^score_')))) |>
    mutate(across(matches('^score_'), function(x) x / temp_sum)) |>
    select(-temp_sum)

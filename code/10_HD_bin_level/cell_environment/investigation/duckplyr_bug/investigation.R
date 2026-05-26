library(tidyverse)
library(data.table)
library(here)

k = 10
cluster_raw_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'ficture_outputs', 'cleaningy', 'k_%d', 'analysis', 'nF%d.d_12',
    'cleaningy_joined_input.tsv.gz'
) |> sprintf(k, k)
out_path = '/fastscratch/myscratch/neagles/smaller.tsv.gz'
max_rows = 1e5

ficture_df = fread(cluster_raw_path, nrows = max_rows, sep = '\t') |>
    as_tibble() |>
    write_tsv(out_path)

a = ficture_df |>
    distinct(barcode, sample_id, factor_K1) |>
    write_tsv('/fastscratch/myscratch/neagles/dplyr.csv')

library(duckplyr)
b = read_csv_duckdb(out_path) |>
    distinct(barcode, sample_id, factor_K1) |>
    collect() |>
    write_tsv('/fastscratch/myscratch/neagles/duckplyr.csv')

# Then reload R
library(tidyverse)
library(here)
library(testthat)

a = read_tsv('/fastscratch/myscratch/neagles/dplyr.csv')
b = read_tsv('/fastscratch/myscratch/neagles/duckplyr.csv')

d = inner_join(
    a |> dplyr::rename(a = factor_K1),
    b |> dplyr::rename(b = factor_K1),
    by = c('barcode', 'sample_id')
)
expect_identical(d$a, d$b)

ref_spe = a |>
    arrange(sample_id, barcode) |>
    slice_sample(prop = 1) |>
    select(sample_id, barcode)

a_clus = ref_spe |>
    left_join(a |> slice_sample(prop = 1), by = c('sample_id', 'barcode')) |>
    pull(factor_K1)

library(duckplyr)
b_clus = ref_spe |>
    left_join(b, by = c('sample_id', 'barcode')) |>
    pull(factor_K1)

expect_identical(a_clus, b_clus)

# After thinking for a while, I'm highly suspicious that maybe in the larger
# real data, left_join() is not preserving row order through duckplyr. In
# this smaller case it appears to though
# https://duckplyr.tidyverse.org/articles/limits.html#output-order-stability

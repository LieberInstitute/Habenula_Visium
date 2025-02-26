#   Before we have actual 2um FICTURE data joined with barcodes, generate a
#   synthetic version so I can write the downstream bin2cell code almost as is
#   in advance

library(here)
library(data.table)
library(tidyverse)
library(sessioninfo)

ficture_input_path = here(
    'processed-data', '10_HD_bin_level', 'ficture', 'inputs',
    'H1-MVPY9BW_A1_8433', 'transcripts.sorted.tsv.gz'
)
positions_path = '/fastscratch/myscratch/neagles/H1-MVPY9BW_A1_8433_tissue_positions.csv.gz'
out_path = here(
    'processed-data', '09_HD_cell_level', 'ficture_aggregate', 'synthetic_input.tsv.gz'
)

dir.create(dirname(out_path), showWarnings = FALSE)

positions = fread(positions_path) |>
    as_tibble()
ficture_input = fread(ficture_input_path) |>
    as_tibble()

#   Randomly assign a barcode to each unique coordinate
temp = ficture_input |>
    group_by(X, Y) |>
    slice_head(n = 1) |>
    select(X, Y) |>
    ungroup()
temp$barcode = positions$barcode[seq_len(nrow(temp))]

left_join(ficture_input, temp, by = c("X", "Y")) |>
    write_tsv(out_path)

session_info()

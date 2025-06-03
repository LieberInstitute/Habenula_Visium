library(here)
library(scran)
library(tidyverse)
library(SpatialExperiment)
library(rjson)
library(sessioninfo)

sample_id_path = here('raw-data', 'sample_info', 'hd_sample_list.txt')
spe_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'cell_environment',
    'spe_filtered.rds'
)
scalefactors_path = here(
    'processed-data', '01_spaceranger', 'probe_fix', '%s', 'outs',
    'binned_outputs', 'square_002um', 'spatial', 'scalefactors_json.json'
)
counts_out_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'cell_environment',
    'ficture_inputs', 'cleaningy_input.tsv.gz'
)
minmax_out_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'cell_environment',
    'ficture_inputs', 'cleaningy_minmax.tsv'
)
out_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'cell_environment',
    'cleaningY', 'temp_chunks', '%d.rds'
)
buffer_prop = 0.05
num_chunks = 50

dir.create(dirname(counts_out_path), showWarnings = FALSE)

################################################################################
#   Read in batch-corrected counts
################################################################################

message(Sys.time(), ' | Reading in SPE and batch-corrected counts...')

spe = readRDS(spe_path)

#   Read in cleaningY chunks of counts and merge
cleaned_counts = list()
for (i in seq_len(num_chunks)) {
    cleaned_counts[[i]] = readRDS(sprintf(out_path, i))
}
cleaned_counts = do.call(rbind, cleaned_counts)

#   Attach to SPE
stopifnot(identical(dim(spe), dim(cleaned_counts)))
assays(spe)$normcounts = cleaned_counts
rm(cleaned_counts)
gc()

################################################################################
#   Convert counts to FICTURE input
################################################################################

sample_ids = readLines(sample_id_path)

#-------------------------------------------------------------------------------
#   Form basic tibble
#-------------------------------------------------------------------------------

#   Form a tibble of the nonzero elements of the normalized-counts matrix
#   (in a memory-efficient way)
message(Sys.time(), ' | Converting counts to FICTURE input...')
counts_mat = as(assays(spe)$normcounts, "TsparseMatrix")
counts_df = tibble(
    X = spatialCoords(spe)[counts_mat@j + 1, 1],
    Y = spatialCoords(spe)[counts_mat@j + 1, 2],
    gene = rownames(spe)[counts_mat@i + 1],
    Count = counts_mat@x,
    sample_id = spe$sample_id[counts_mat@j + 1],
    barcode = colnames(spe)[counts_mat@j + 1]
)

#   Convert units of spatial coords to microns
for (sample_id in sample_ids) {
    micron_per_px = fromJSON(file = sprintf(scalefactors_path, sample_id))[['microns_per_pixel']]
    counts_df[counts_df$sample_id == sample_id, 'X'] = counts_df[counts_df$sample_id == sample_id, 'X'] * micron_per_px
    counts_df[counts_df$sample_id == sample_id, 'Y'] = counts_df[counts_df$sample_id == sample_id, 'Y'] * micron_per_px
}

#-------------------------------------------------------------------------------
#   Put each samples next to each other
#-------------------------------------------------------------------------------

#   Find a range of X values slightly larger than any particular sample
x_size = counts_df |>
    group_by(sample_id) |>
    summarize(x_diff = max(X) - min(X)) |>
    pull(x_diff) |>
    max()
x_size = (1 + buffer_prop) * x_size

#   Separate samples by placing each sample into the same Y range and adjacent
#   X ranges
counts_df = counts_df |>
    group_by(sample_id) |>
    mutate(
        X = X - min(X) + x_size * (match(cur_group()$sample_id, sample_ids) - 1),
        Y = Y - min(Y)
    ) |>
    ungroup() |>
    arrange(X)

#-------------------------------------------------------------------------------
#   Write FICTURE inputs
#-------------------------------------------------------------------------------

write_tsv(counts_df, counts_out_path)

#   Also write coordinate ranges to a TSV file (another required FICTURE input)
counts_df |>
    summarize(xmin = min(X), xmax = max(X), ymin = min(Y), ymax = max(Y)) |>
    t() |>
    as.data.frame() |>
    rownames_to_column('coord_type') |>
    write_tsv(minmax_out_path, col_names = FALSE)

session_info()

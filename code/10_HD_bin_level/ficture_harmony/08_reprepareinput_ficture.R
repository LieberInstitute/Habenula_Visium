library(sessioninfo)
library(here)
library(SpatialExperiment)
library(scran)
library(scater)
library(tidyverse)
library(rjson)

sample_id_path = here('raw-data', 'sample_info', 'hd_sample_list.txt')
spe_in_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'ficture_harmony','spe',
    'y_clean_spe.rds'
)
scalefactors_path = here(
    'processed-data', '01_spaceranger', 'probe_fix', '%s', 'outs',
    'binned_outputs', 'square_002um', 'spatial', 'scalefactors_json.json'
)
counts_out_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'ficture_harmony',
    'ficture_inputs', 'cleany', 'input.tsv.gz'
)
minmax_out_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'ficture_harmony',
    'ficture_inputs', 'cleany', 'minmax.tsv'
)
buffer_prop = 0.05

dir.create(dirname(counts_out_path), showWarnings = FALSE)

spe = readRDS(spe_in_path)
sample_ids = readLines(sample_id_path)

#   Filter raw SPE: drop bins with 0 counts for all genes, and drop genes with
#   0 counts in every bin
message(Sys.time(), ' | Filtering bins and genes...')
spe = spe[rowSums(assays(spe)$counts) > 0, colSums(assays(spe)$counts) > 0]

#   Form a tibble of the nonzero elements of the normalized-counts matrix
#   (in a memory-efficient way)
message(Sys.time(), ' | Converting counts to FICTURE input...')
counts_mat = as(assays(spe)$counts, "TsparseMatrix")
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

write_tsv(counts_df, counts_out_path)

#   Also write coordinate ranges to a TSV file (another required FICTURE input)
counts_df |>
    summarize(xmin = min(X), xmax = max(X), ymin = min(Y), ymax = max(Y)) |>
    t() |>
    as.data.frame() |>
    rownames_to_column('coord_type') |>
    write_tsv(minmax_out_path, col_names = FALSE)

session_info()
print("1")

library(here)
library(tidyverse)
library(SpatialExperiment)
library(HDF5Array)
library(crawdad)
library(rjson)
library(sessioninfo)

spe_dir = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'spe_norm_filtered'
)
banksy_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'banksy', 'lambda0_2',
    'leiden_res1_3_subset.csv'
)
scalefactors_path = here(
    'processed-data', '01_spaceranger', 'probe_fix', '%s', 'outs',
    'binned_outputs', 'square_002um', 'spatial', 'scalefactors_json.json'
)
sample_id_path = here('raw-data', 'sample_info', 'hd_sample_list.txt')
plot_path = here(
    'plots', '09_HD_cell_level', 'probe_fix', 'crawdad', 'dot_plot.pdf'
)
# scales = c(100, 200, 500, 1000, 5000)
scales = c(500, 5000)
random_seed = 0
downsample_prop = 0.1

dir.create(dirname(plot_path), showWarnings = FALSE)

sample_id = readLines(sample_id_path)[
    as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
]

spe = loadHDF5SummarizedExperiment(spe_dir)
spe = spe[, spe$sample_id == sample_id]

#   Ultimately, we'll be converting spatial coordinates to units of microns,
#   which is more interpretable than pixels
micron_per_px = fromJSON(file = sprintf(scalefactors_path, sample_id))[['microns_per_pixel']]

#   Gather spatial coordinates and Banksy clusters
cell_df = tibble(
        key = spe$key,
        x = spatialCoords(spe)[, 'pxl_col_in_fullres'] * micron_per_px,
        y = spatialCoords(spe)[, 'pxl_row_in_fullres'] * micron_per_px
    ) |>
    left_join(read_csv(banksy_path, show_col_types = FALSE), by = 'key') |>
    slice_sample(prop = downsample_prop) |>
    mutate(banksy = factor(banksy_lambda0_2)) |>
    select(x, y, banksy) |>
    as.data.frame()
stopifnot(!any(is.na(cell_df$banksy)))

pos_df = toSF(pos = select(cell_df, c(x, y)), cellTypes = cell_df$banksy)

#   Shuffle cell-type assignments to create null background
shuffle_list = makeShuffledCells(
    pos_df, scales = scales, seed = random_seed, verbose = TRUE
)

#   Main Z-score calculation for each reference-neighbor pair
results = findTrends(
    pos_df, shuffleList = shuffle_list, verbose = TRUE, returnMeans = FALSE
)

#   Reformat and compute multiple-testing-corrected Z-score
results = meltResultsList(results, withPerms = TRUE)
z_sig = correctZBonferroni(results)

#   Main dot plot figure
pdf(plot_path)
vizColocDotplot(
    results, zSigThresh = z_sig, zScoreLimit = 2 * z_sig, dotSizes = c(1, 5)
)
dev.off()

session_info()

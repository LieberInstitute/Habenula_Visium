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


sample_id = readLines(sample_id_path)[
    as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
]

spe = loadHDF5SummarizedExperiment(spe_dir)
spe = spe[, spe$sample_id == sample_id]

#   Join in Banksy clusters with SPE
spe$banksy = tibble(key = spe$key) |>
    left_join(read_csv(banksy_path, show_col_types = FALSE), by = 'key') |>
    pull(banksy_lambda0_2)
stopifnot(!any(is.na(spe$banksy)))

#   Get spatial coordinates in units of microns
micron_per_px = fromJSON(file = sprintf(scalefactors_path, sample_id))[['microns_per_pixel']]
coords_df = spatialCoords(spe) |>
    as.data.frame() |>
    mutate(
        x = pxl_col_in_fullres * micron_per_px,
        y = pxl_row_in_fullres * micron_per_px
    ) |>
    select(x, y)

pos_df = toSF(pos = coords_df, cellTypes = factor(spe$banksy))

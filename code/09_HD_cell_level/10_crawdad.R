library(here)
library(tidyverse)
library(SpatialExperiment)
library(HDF5Array)
library(crawdad)
library(sessioninfo)

spe_dir = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'spe_norm_filtered'
)
banksy_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'banksy', 'lambda0_2',
    'leiden_res1_3_subset.csv'
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

coords_df = as.data.frame(spatialCoords(spe))
colnames(coords_df) = c('x', 'y')

pos_df = toSF(pos = coords_df, cellTypes = factor(spe$banksy))

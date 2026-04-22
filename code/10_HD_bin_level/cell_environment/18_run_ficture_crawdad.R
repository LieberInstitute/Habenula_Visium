library(here)
library(tidyverse)
library(crawdad)
library(sessioninfo)

sample_info_path = here('raw-data', 'sample_info', 'hd_basic_info_split.csv')
sample_info = read_csv(sample_info_path, show_col_types = FALSE)
this_sample_id = sample_info$tissue_id[
    as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
]

in_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'crawdad', 'input_cells.csv.gz'
)
out_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'crawdad', 'output', sprintf('%s_results.csv', this_sample_id)
)
scales = c(100, 200, 500, 1000, 5000)
random_seed = 0

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
dir.create(dirname(out_path), showWarnings = FALSE)

cell_df = read_csv(in_path, show_col_types = FALSE) |>
    filter(sample_id == this_sample_id) |>
    #   Pretend here that we have 2 cells: one for the Banksy cluster and one
    #   for the FICTURE (same coordinates)
    mutate(ficture_cluster = paste0('Factor_', ficture_cluster)) |>
    pivot_longer(
        c(cell_type, ficture_cluster), names_to = NULL, values_to = 'cluster'
    ) |>
    as.data.frame()

pos_df = toSF(pos = select(cell_df, c(x, y)), cellTypes = cell_df$cluster)

#   Shuffle cell-type assignments to create null background
shuffle_list = makeShuffledCells(
    pos_df, scales = scales, seed = random_seed, ncores = num_cores,
    verbose = TRUE
)

#   Main Z-score calculation for each reference-neighbor pair
results = findTrends(
        pos_df, shuffleList = shuffle_list, returnMeans = FALSE,
        ncores = num_cores, verbose = TRUE
    ) |>
    #   Reformat and export
    meltResultsList(withPerms = TRUE) |>
    write_csv(out_path)

message("Memory usage:")
gc()

session_info()

## This script was made using slurmjobs version 1.3.0
## available from http://research.libd.org/slurmjobs/

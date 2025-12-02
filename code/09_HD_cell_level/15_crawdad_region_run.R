library(getopt)
library(here)
library(tidyverse)
library(crawdad)
library(sessioninfo)

# Import command-line parameters
spec <- matrix(
    c(
        c("sample_id", "region"),
        c("s", "r"),
        rep("1", 2),
        rep("character", 2),
        rep("Add variable description here", 2)
    ),
    ncol = 5
)
opt <- getopt(spec)

message("Using the following parameters:")
print(opt)

in_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'crawdad', 'region',
    'input_cells.csv.gz'
)
out_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'crawdad', 'region',
    'output', sprintf('%s_%s_results.csv', opt$sample_id, opt$region)
)
scales = c(100, 200, 500, 1000, 5000)
random_seed = 0

num_cores = as.integer(Sys.getenv("SLURM_CPUS_ON_NODE"))
dir.create(dirname(out_path), showWarnings = FALSE)

cell_df = read_csv(in_path, show_col_types = FALSE) |>
    filter(sample_id == opt$sample_id, region_anno == opt$region, !drop) |>
    as.data.frame()

pos_df = toSF(pos = select(cell_df, c(x, y)), cellTypes = cell_df$cell_type)

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

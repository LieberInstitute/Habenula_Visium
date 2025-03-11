library(here)
library(tidyverse)
library(scran)
library(spatialLIBD)
library(sessioninfo)
library(HDF5Array)

sample_id_path = here('raw-data', 'sample_info', 'hd_sample_list.txt')
spe_out_path = here(
    'processed-data', '10_HD_bin_level', 'ficture_harmony', 'spe_raw.rds'
)

sample_ids = readLines(sample_id_path)
sr_out_dirs = here(
    'processed-data', '01_spaceranger', sample_ids, 'outs',
    'binned_outputs', 'square_002um'
)

#   Hack around 'read10xVisium's requirement for the 'outs' directory to be the
#   immediate parent to 'spatial' directory and other outputs (create a symlink
#   named 'outs'). Note we already handled the other required workaround: create
#   a tissue_positions.csv file (not parquet format)
temp_sr_dirs = file.path(tempdir(), sample_ids, 'outs')
for (this_dir in temp_sr_dirs) {
    dir.create(dirname(this_dir))
}

file.symlink(sr_out_dirs, temp_sr_dirs) |>
    all() |>
    stopifnot()

message(Sys.time(), ' | Building SpatialExperiment...')
spe <- read10xVisium(
    samples = temp_sr_dirs,
    sample_id = sample_ids,
    type = "sparse",
    data = "raw",
    images = "lowres",
    load = FALSE
)

message(Sys.time(), " | Saving raw SPE")
saveRDS(spe, spe_out_path)

session_info()

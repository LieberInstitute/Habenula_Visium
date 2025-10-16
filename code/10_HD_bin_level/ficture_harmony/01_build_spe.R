#   Build a basic 2um QC'd bin-level SpatialExperiment object, which downstream
#   will eventually become used as input for FICTURE

library(here)
library(tidyverse)
library(scran)
library(spatialLIBD)
library(sessioninfo)
library(HDF5Array)

sample_info_path = here('raw-data', 'sample_info', 'hd_basic_info.csv')
sample_info = read_csv(sample_info_path, show_col_types = FALSE)

spe_out_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples', 'ficture_harmony',
    'spe_raw.rds'
)
sr_out_dirs = here(
    sample_info$spaceranger_dir, 'outs', 'binned_outputs', 'square_002um'
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
    load = TRUE
)

#   Drop problematic tissue regions found in the 8um bin-level QC, other than
#   one UMI-based filter (the used cutoff of 5 at 8um is so low that there is no
#   good equivalent at 2um; FICTURE should drop this region anyway). In an
#   interactive test, I showed that the same tissue regions are dropped at 2um
#   when array coordinates are multiplied by 4 relative to 8um
message(Sys.time(), ' | Filtering bins according to 8um bin-level QC...')
spe = spe[, (spe$sample_id != 'H1-MVPY9BW_A1_8433') | (spe$array_col <= 793 * 4)]
spe = spe[, (spe$sample_id != 'H1-XQQD7C7_A1_8518') | (spe$array_row <= 764 * 4)]
spe = spe[, (spe$sample_id != 'H1-XQQD7C7_A1_8518') | (spe$array_col > 142 * 4)]

message(Sys.time(), " | Saving raw SPE")
saveRDS(spe, spe_out_path)

session_info()

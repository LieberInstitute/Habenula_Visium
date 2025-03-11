library(sessioninfo)
library(here)
library(HDF5Array)
library(SpatialExperiment)
library(harmony)

spe_in_dir = here(
    'processed-data', '10_HD_bin_level', 'ficture_harmony', 'spe_raw'
)

set.seed(0)

#   Load and run Harmony
spe = loadHDF5SummarizedExperiment(spe_in_dir)
spe = RunHarmony(spe, group.by.vars = "sample_id", dims.use = "PCA")

#   Save embedding in place
quickResaveHDF5SummarizedExperiment(spe)

session_info()

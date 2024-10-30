library(here)
library(tidyverse)
library(spatialLIBD)
library(sessioninfo)
library(HDF5Array)
library(SEraster)
library(rjson)

#   Number of times more bins than a Visium experiment (per dimension)
res_scalar = 1

sample_id = 'H1-W369TJK_D1_9090'
spe_norm_dir = here('processed-data', '10_HD_bin_level', 'spe_norm')
spe_out_dir = here(
    'processed-data', '10_HD_bin_level', 'rasterized',
    sprintf('spe_%s_%sx_standard_res', sample_id, res_scalar)
)
json_path = here(
    'processed-data', '01_spaceranger', sample_id, 'outs', 'binned_outputs',
    'square_008um', 'spatial', 'scalefactors_json.json'
)

num_cores = as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))
dir.create(dirname(spe_out_dir), showWarnings = FALSE)

#   Calculate the appropriate SEraster resolution (measured in pixels). The
#   "100" here represents the 100um inter-spot distance in Visium standard
res = round(
    100 / fromJSON(file = json_path)[['microns_per_pixel']] / res_scalar
)

#   Load and subset to this sample
message(Sys.time(), ' | Loading and bringing into memory')
spe = loadHDF5SummarizedExperiment(spe_norm_dir)

#   Bring into memory to speed up computations
assays(spe)$logcounts = as(assays(spe)$logcounts, "dgCMatrix")

#   Rasterize
message(
    Sys.time(),
    sprintf(
        " | Rasterizing with %s-pixel width, which is %s times the resolution (per dimension) of Visium standard.",
        res, res_scalar
    )
)
spe_raster = rasterizeGeneExpression(
    spe, assay_name = "logcounts", resolution = res, n_threads = num_cores
)

#   Add back rowData, which is removed for some reason
stopifnot(identical(rownames(spe_raster), rownames(spe)))
rowData(spe_raster) = rowData(spe)

#   Save
message(Sys.time(), ' | Saving rasterized SPE')
spe_raster <- saveHDF5SummarizedExperiment(
    spe_raster, dir = spe_out_dir, replace = TRUE, as.sparse = TRUE
)

session_info()

## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/

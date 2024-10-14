library(here)
library(tidyverse)
library(scran)
library(spatialLIBD)
library(sessioninfo)
library(HDF5Array)

sample_id = 'H1-W369TJK_D1_9090'
spe_raw_dir = here('processed-data', '10_HD_bin_level', 'spe_raw')
spe_norm_dir = here('processed-data', '10_HD_bin_level', 'spe_norm')
plot_dir = here('plots', '06_bin_level')
wm_genes = c("MBP", "GFAP", "PLP1", "AQP4")
sr_out_dir = here(
    'processed-data', '01_spaceranger', sample_id, 'outs',
    'binned_outputs', 'square_008um'
)
reference_gtf = '/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf'

#   Hack around 'read10xVisium's requirement for the 'outs' directory to be the
#   immediate parent to 'spatial' directory and other outputs (create a symlink
#   named 'outs'). Note we already handled the other required workaround: create
#   a tissue_positions.csv file (not parquet format)
temp_sr_dir = file.path(tempdir(), sample_id, 'outs')
dir.create(dirname(temp_sr_dir))

file.symlink(sr_out_dir, temp_sr_dir) |>
    all() |>
    stopifnot()

#   Note providing the reference GTF here is mandatory (without doing so,
#   read10xVisiumWrapper searches for a web summary that doesn't exist to try
#   to infer the GTF)
message(Sys.time(), ' | Building SpatialExperiment...')
spe <- read10xVisiumWrapper(
    samples = temp_sr_dir,
    sample_id = sample_id,
    type = "sparse",
    data = "raw",
    images = "lowres",
    load = FALSE,
    reference_gtf = reference_gtf
)

message(Sys.time(), " - Saving raw SPE")
spe <- saveHDF5SummarizedExperiment(
    spe, dir = spe_raw_dir, replace = TRUE, as.sparse = TRUE
)

#   Filter raw SPE: take only bins in tissue, drop bins with 0 counts for all
#   genes, and drop genes with 0 counts in every bin
message(Sys.time(), ' | Filtering bins, genes, and to tissue...')
spe <- spe[
    rowSums(assays(spe)$counts) > 0,
    (colSums(assays(spe)$counts) > 0) & spe$in_tissue
]

#   Use library-size normalization (normalization by deconvolution is not
#   computationally feasible with data this large)
message(Sys.time(), ' | Performing log normalization...')
spe = computeLibraryFactors(spe)
spe = logNormCounts(spe)

message(Sys.time(), " - Saving normalized SPE")
spe <- saveHDF5SummarizedExperiment(
    spe, dir = spe_norm_dir, replace = TRUE, as.sparse = TRUE
)

session_info()

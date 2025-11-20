#   Build a 8um bin-level SpatialExperiment object with basic expression
#   filtering and library-size normalization

library(here)
library(tidyverse)
library(scran)
library(spatialLIBD)
library(sessioninfo)
library(HDF5Array)

sample_info_path = here('raw-data', 'sample_info', 'hd_basic_info.csv')
spe_raw_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'spe_raw.rds'
)
spe_raw_old_dir = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'spe_raw'
)
spe_norm_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'spe_norm.rds'
)

sample_info = read_csv(sample_info_path, show_col_types = FALSE)
sample_ids = sample_info$sample_id
sr_out_dirs = here(
    sample_info$spaceranger_dir, 'outs', 'binned_outputs', 'square_008um'
)
reference_gtf = '/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz'

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

#   Note providing the reference GTF here is mandatory (without doing so,
#   read10xVisiumWrapper searches for a web summary that doesn't exist to try
#   to infer the GTF)
message(Sys.time(), ' | Building SpatialExperiment...')
spe <- read10xVisiumWrapper(
    samples = temp_sr_dirs,
    sample_id = sample_ids,
    type = "sparse",
    data = "raw",
    images = "lowres",
    load = TRUE,
    reference_gtf = reference_gtf
)

message(Sys.time(), " | Saving raw SPE")
saveRDS(spe, spe_raw_path)

#   For some reason, 'in_tissue' has incorrect values in the latest spaceranger
#   run, with some truly in-tissue regions marked FALSE. Only the samples from
#   batch 1 and 2 were badly problematic, so we just take 'in_tissue' from the
#   old SPE for those samples
spe_old = loadHDF5SummarizedExperiment(spe_raw_old_dir)
spe_old = spe_old[, spe_old$sample_id %in% sample_ids]
stopifnot(all(spe_old$key %in% spe$key))
spe$in_tissue[match(spe_old$key, spe$key)] = spe_old$in_tissue

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

message(Sys.time(), " | Saving normalized SPE")
saveRDS(spe, spe_norm_path)

session_info()

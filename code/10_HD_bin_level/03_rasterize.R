library(here)
library(tidyverse)
library(spatialLIBD)
library(sessioninfo)
library(SEraster)
library(rjson)
library(scran)

#   Number of times more bins than a Visium experiment (per dimension)
res_scalar = 1

sample_info_path = here('raw-data', 'sample_info', 'hd_basic_info.csv')
sample_info = read_csv(sample_info_path)
task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
sample_id = sprintf(
    'Br%s_%d',
    str_extract(sample_info$sample_id[(task_id + 1) %/% 2], '[0-9]{4}$'),
    (task_id - 1) %% 2 + 1
)
spaceranger_dir = sample_info$spaceranger_dir[(task_id + 1) %/% 2]

spe_norm_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
spe_out_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'rasterized',
    sprintf('spe_%s_lowres.rds', sample_id)
)
json_path = here(
    spaceranger_dir, 'outs', 'binned_outputs', 'square_008um', 'spatial',
    'scalefactors_json.json'
)
plot_path = here(
    'plots', '10_HD_bin_level', 'no_secondary', 'rasterized',
    sprintf('WM_%s.pdf', sample_id)
)
markers = c("MBP", "GFAP", "PLP1", "AQP4")

num_cores = as.integer(Sys.getenv("SLURM_CPUS_ON_NODE"))
dir.create(dirname(spe_out_path), showWarnings = FALSE)
dir.create(dirname(plot_path), showWarnings = FALSE)

#   Calculate the appropriate SEraster resolution (measured in pixels). The
#   "100" here represents the 100um inter-spot distance in Visium standard
res = round(
    100 / fromJSON(file = json_path)[['microns_per_pixel']] / res_scalar
)

#   Load and subset to this sample
message(Sys.time(), ' | Loading this sample and bringing into memory')
spe = readRDS(spe_norm_path)
spe = spe[, spe$sample_id == sample_id]

#   Rasterize
message(
    Sys.time(),
    sprintf(
        " | Rasterizing with %s-pixel width, which is %s times the resolution (per dimension) of Visium standard.",
        res, res_scalar
    )
)
spe_raster = rasterizeGeneExpression(
    spe, assay_name = "counts", resolution = res, n_threads = num_cores,
    fun = "sum"
)

#   Fix several object attributes
names(assays(spe_raster)) = 'counts'
colnames(spatialCoords(spe_raster)) = c('pxl_col_in_fullres', 'pxl_row_in_fullres')
spe_raster$sample_id = sample_id

#   Add back imgData and rowData, which are removed for some reason
stopifnot(identical(rownames(spe_raster), rownames(spe)))
rowData(spe_raster) = rowData(spe)
imgData(spe_raster) = imgData(spe)

#   Log-normalize rasterized raw counts
message(Sys.time(), ' | Performing log normalization...')
spe_raster = computeLibraryFactors(spe_raster)
spe_raster = logNormCounts(spe_raster)

#   To visually verify we're getting reasonable results, plot some markers on
#   the rasterized data
markers = rownames(spe_raster)[match(markers, rowData(spe_raster)$gene_name)]
pdf(plot_path)
vis_gene(
    spe_raster, geneid = markers, multi_gene_method = "z_score",
    auto_crop = FALSE, point_size = 1, spatial = TRUE
)
dev.off()

#   Save
message(Sys.time(), ' | Saving rasterized SPE')
saveRDS(spe_raster, spe_out_path)

session_info()

## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/

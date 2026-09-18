library(SpatialExperiment)
library(here)
library(tidyverse)
library(sessioninfo)
library(Matrix)
library(nnSVG)

sample_info_path = here('raw-data', 'sample_info', 'hd_basic_info.csv')
sample_info = read_csv(sample_info_path, show_col_types = FALSE)
task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
sample_id = sprintf(
    'Br%s_%d',
    str_extract(sample_info$sample_id[(task_id + 1) %/% 2], '[0-9]{4}$'),
    (task_id - 1) %% 2 + 1
)

spe_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'rasterized',
    sprintf('spe_%s_lowres.rds', sample_id)
)
out_path <- here(
    "processed-data", '10_HD_bin_level', 'no_secondary', "nnSVG_out",
    paste0(sample_id, ".csv")
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_ON_NODE"))
set.seed(0)
dir.create(dirname(out_path), showWarnings = FALSE)

#-------------------------------------------------------------------------------
#   Subset to this sample
#-------------------------------------------------------------------------------

message(Sys.time(), " | Loading, subsetting, and bringing assays into memory")
spe <- readRDS(spe_path)
spe <- spe[, spe$sample_id == sample_id]

#-------------------------------------------------------------------------------
#   Filter lowly expressed and mitochondrial genes, and low-count spots
#-------------------------------------------------------------------------------

message(Sys.time(), " | Filtering genes and spots")
spe <- filter_genes(
    spe,
    filter_genes_ncounts = 3,
    filter_genes_pcspots = 5,
    filter_mito = TRUE
)

#   This step is not in the vignette but seems to eliminate BRISC estimation
#   errors (see https://github.com/lmweber/nnSVG/issues/16)
enough_counts = colSums(assays(spe)$counts) >= 10
message(
    sprintf(
        "Dropping %.1f%% of spots with low counts", 
        100 * (1 - mean(enough_counts))
    )
)
spe = spe[, enough_counts]

message("Dimensions of spe after filtering:")
print(dim(spe))

#-------------------------------------------------------------------------------
#   Run nnSVG and export results
#-------------------------------------------------------------------------------

message(Sys.time(), " | Running nnSVG")
spe <- nnSVG(spe, n_threads = num_cores)

message(Sys.time(), " | Exporting results")
write_csv(as_tibble(rowData(spe)), out_path)

session_info()

library(tidyverse)
library(here)
library(sessioninfo)
library(SingleCellExperiment)
library(spatialLIBD)
library(edgeR)
library(limma)
library(duckplyr)

#   Grab cell type from array task ID
task_id = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
sce_dir = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'cell_vs_extra_DE',
    'sce'
)
cell_type_clean = list.files(sce_dir, pattern = '^sce_[^a].*_pb\\.rds$')[task_id] |> 
    str_extract('^sce_(.*)_pb.rds$', group = 1)
cell_type = cell_type_clean |> 
    str_replace_all('--', '/') |> 
    str_replace_all('-', '.')

sce_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'cell_vs_extra_DE',
    'sce', sprintf('sce_%s_pb.rds', cell_type_clean)
)
out_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'cell_vs_extra_DE',
    'main_results', sprintf('DE_%s.parquet', cell_type_clean)
)
covariates = c('compartment', 'sample_id')

dir.create(dirname(out_path), showWarnings = FALSE)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

sce = readRDS(sce_path)
sce$sample_id = factor(sce$sample_id, levels = sort(unique(sce$sample_id)))
sce$donor = factor(sce$donor, levels = sort(unique(sce$donor)))

des = paste('~', paste0(covariates, collapse = ' + ')) |>
    as.formula() |>
    model.matrix(data = colData(sce)) |>
    as.data.frame()

dge = calcNormFactors(sce)
dge = dge[filterByExpr.DGEList(dge, design = des), , keep.lib.sizes = FALSE]
dge = calcNormFactors(dge)

#   DE
de_df = dge |>
    voomLmFit(
        design = des, adaptive.span = TRUE, sample.weights = TRUE,
        block = dge$donor
    ) |>
    eBayes() |>
    topTable(coef = "compartmentextra", number = Inf) |>
    as_tibble() |>
    mutate(cell_type = cell_type) |>
    select(
        gene_id, gene_name, cell_type, logFC, AveExpr, t, P.Value, adj.P.Val, B
    ) |>
    compute_parquet(out_path)

session_info()

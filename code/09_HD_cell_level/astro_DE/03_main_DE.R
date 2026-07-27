library(tidyverse)
library(here)
library(sessioninfo)
library(qs2)
library(SingleCellExperiment)
library(spatialLIBD)
library(edgeR)
library(limma)
library(duckplyr)

task_id = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

sce_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'astro_DE',
    'astro_sce.qs2'
)
out_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'astro_DE',
    'main_results', sprintf('DE_%d.parquet', task_id)
)
dge_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'astro_DE',
    'main_results', 'DGE.qs2'
)
cont_covariates = c('ncells', 'expr_chrM_ratio')

set.seed(task_id)
dir.create(dirname(out_path), showWarnings = FALSE)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

sce = qs_read(sce_path)

#   task_id = 0 is reserved for the true DE of the proper labels. All other
#   tasks permute labels to as a whole generate a null distribution of t-stats
if (task_id > 0) {
    message("Permuting labels...")
    sce$astro_label = tibble(astro_label = sce$astro_label) |>
        slice_sample(n = ncol(sce)) |>
        pull(astro_label)
} else {
    message("Using original labels.")
}
sce = sce[, sce$astro_label %in% c('medial', 'lateral')]

#   Pseudobulk by astro_label and tissue ID
sce_pb = registration_pseudobulk(
    sce, var_registration = "astro_label", var_sample_id = "sample_id"
)

#   Recompute important colData and clean things up
sce_pb$pb_sample_id = colnames(sce_pb)
sce_pb$donor = sub('_[12]$', '', sce_pb$sample_id)
sce_pb$sum_umi = unname(colSums(counts(sce_pb)))
sce_pb$expr_chrM = colSums(
    counts(sce_pb)[which(seqnames(sce_pb) == "chrM"), , drop = FALSE]
)
sce_pb$expr_chrM_ratio = sce_pb$expr_chrM / sce_pb$sum_umi
sce_pb$astro_label = factor(sce_pb$astro_label, levels = c('medial', 'lateral'))
sce_pb$donor = factor(sce_pb$donor, levels = sort(unique(sce_pb$donor)))

#   Center and scale continuous covariates
for (this_covariate in cont_covariates) {
    sce_pb[[this_covariate]] = as.numeric(scale(sce_pb[[this_covariate]]))
}

des = paste('~ astro_label +', paste0(cont_covariates, collapse = ' + ')) |>
    as.formula() |>
    model.matrix(data = colData(sce_pb)) |>
    as.data.frame()

dge = calcNormFactors(sce_pb)
dge = dge[filterByExpr.DGEList(dge, design = des), , keep.lib.sizes = FALSE]
dge = calcNormFactors(dge)

#   DE
dge = dge |>
    voomLmFit(
        design = des, adaptive.span = TRUE, sample.weights = TRUE,
        block = dge$donor
    )

de_df = dge |>
    eBayes() |>
    topTable(coef = "astro_labellateral", number = Inf) |>
    as_tibble()

#   For the real DE, save additional info. Otherwise, we really only need gene
#   and t-stats for the permutations
if (task_id == 0) {
    de_df |>
        select(gene_id, gene_name, t, logFC, adj.P.Val) |>
        compute_parquet(out_path)

    qs_save(dge, dge_path)
} else {
    de_df |>
        select(gene_id, t) |>
        mutate(permutation = task_id) |>
        compute_parquet(out_path)
}

session_info()

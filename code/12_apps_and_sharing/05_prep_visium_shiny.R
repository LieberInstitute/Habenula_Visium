library(here)
library(tidyverse)
library(sessioninfo)
library(SpatialExperiment)
library(qs2)
library(spatialLIBD)
library(lobstr)
library(withr)

spe_in_path = here("processed-data", "04_harmony_BayesSpace", "spe_harmony.rds")
sce_pb_in_path = here(
    "processed-data", "05_brain_area_differential_expression",
    "sce_pseudo_BayesSpace_k09.rds"
)
cluster_dir = here(
    "processed-data", "04_harmony_BayesSpace", "clusters_BayesSpace"
)
out_dir = here("processed-data", "12_apps_and_sharing", "05_prep_visium_shiny")
app_dir = here("code", "12_apps_and_sharing", "shiny_visium_app")
spe_col_data_keep = c(
    "sample_id", "key", "array_row", "array_col", "X10x_.*", "sum_umi",
    "sum_gene", "expr_chrM", "expr_chrM_ratio", "ManualAnnotation",
    "brain_id", "brain_area", "age", "sex", "ethnicity", "pmi",
    "rin", "edge_spot", "edge_distance", "SNN_k10_.*",
    "BayesSpace_harmony_.*"
)
sig_genes_n = 1000

dir.create(out_dir, showWarnings = FALSE)
dir.create(app_dir, showWarnings = FALSE)

################################################################################
#   Prep 'spe'
################################################################################

spe = readRDS(spe_in_path)

message("Original 'spe' size:")
print(obj_size(spe))

#   Don't need this for sharing or Shiny
assays(spe)$binomial_deviance_residuals = NULL

spe = cluster_import(spe, cluster_dir = cluster_dir, prefix = "")
colnames(spe) = spe$key

#   There are too many reducedDims, most of which are not needed
reducedDims(spe) = reducedDims(spe)[
    c("PCA", "HARMONY", "UMAP.HARMONY")
]
reducedDimNames(spe) = c("PCA", "PCA_HARMONY", "UMAP_HARMONY")

#   Essentially dropping colData columns that were used for temporary
#   tasks or that contain a single value
col_data_pattern = paste0(
    "^", paste0(spe_col_data_keep, collapse = "$|^"), "$"
)
colData(spe) = colData(spe)[
    , grepl(col_data_pattern, colnames(colData(spe)))
]

################################################################################
#   Prep 'sce_pb'
################################################################################

sce_pb = readRDS(sce_pb_in_path)

sce_pb$BayesSpace_harmony_k09 = sce_pb$BayesSpace
sce_pb$BayesSpace = NULL

colnames(sce_pb) = paste(
    sce_pb$sample_id, sce_pb$BayesSpace_harmony_k09, sep = "_"
)

################################################################################
#   Prep 'sig_genes' for the Shiny app
################################################################################

sce_pb$spatialLIBD = sce_pb$BayesSpace_harmony_k09
sig_genes = sig_genes_extract_all(
    n = min(sig_genes_n, nrow(sce_pb)),
    modeling_results = readRDS(modeling_path), sce_layer = sce_pb
)
sce_pb$spatialLIBD = NULL

################################################################################
#   Export and link into the destination directory
################################################################################

#   For ExperimentHub/ spatialLIBD::fetch_data()
saveRDS(spe, file.path(out_dir, 'spe_visium_habenula_atlas.rds'))
saveRDS(sce_pb, file.path(out_dir, 'sce_pb_visium_habenula_atlas.rds'))

#   For the Shiny app
assays(spe) = list(logcounts = logcounts(spe))
assays(sce_pb) = list(logcounts = logcounts(sce_pb))

qs_save(spe, file.path(out_dir, 'spe_shiny.qs2'))
qs_save(sce_pb, file.path(out_dir, 'sce_pb_shiny.qs2'))
qs_save(sig_genes, file.path(out_dir, 'sig_genes_shiny.qs2'))

for (f_base_name in c('spe_shiny.qs2', 'sce_pb_shiny.qs2', 'sig_genes_shiny.qs2')) {
    with_dir(
    app_dir,
        system(
            sprintf(
                "ln -s ../../processed-data/12_apps_and_sharing/05_prep_visium_shiny/%s %s",
                f_base_name, f_base_name
            )
        )
    )
}

session_info()

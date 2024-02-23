library("spatialLIBD")

spe <- readRDS("spe.rds")

## Quickly explore the data
vars <- colnames(colData(spe))
spatialLIBD::run_app(
    spe = spe,
    sce_layer = NULL,
    modeling_results = NULL,
    sig_genes = NULL,
    spe_discrete_vars = c("ManualAnnotation", "overlaps_tissue", vars[grep("^10x_", vars)]),
    spe_continuous_vars = c("sum_umi", "sum_gene",
        "expr_chrM", "expr_chrM_ratio"),
    default_cluster = "10x_graphclust"
)

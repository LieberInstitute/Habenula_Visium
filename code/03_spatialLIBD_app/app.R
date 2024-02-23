library("spatialLIBD")
library("markdown") ## Hm... to avoid this error
# 2021-11-11T05:30:50.218127+00:00 shinyapps[5096402]: Warning: Error in loadNamespace: there is no package called ‘markdown’

## spatialLIBD uses golem
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

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

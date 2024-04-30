library("spatialLIBD")
library("markdown") ## Hm... to avoid this error
# 2021-11-11T05:30:50.218127+00:00 shinyapps[5096402]: Warning: Error in loadNamespace: there is no package called ‘markdown’

## spatialLIBD uses golem
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())


## I added the symbolic link in the root directory. Fails to read it in the current directory
spe <- readRDS("spe_qc_low_lib_edge.rds") # Note: poner el objeto con los QCs, spe_qc_low_lib_edge.rds / 

# ## Import BayesSpace clusters
# spe <- cluster_import(spe,
#     cluster_dir = "clusters_BayesSpace",
#     prefix = ""
# )

## Quickly explore the data
vars <- colnames(colData(spe))
colnames(colData(spe)) <- vars <- gsub("X10x", "10x", vars)
spatialLIBD::run_app(
    spe = spe,
    sce_layer = NULL,
    modeling_results = NULL,
    sig_genes = NULL,
    spe_discrete_vars = c(
        "ManualAnnotation",
        "overlaps_tissue",
        vars[grep("^10x_", vars)],
        vars[grep("^scran_", vars)],
        "edge_spots"
        #vars[grep("^SNN_k10", vars)],
        #vars[grep("^BayesSpace_harmony_", vars)]
    ),
    spe_continuous_vars = c(
        "sum_umi",
        "sum_gene",
        "expr_chrM",
        "expr_chrM_ratio",
        "edge_distance"
    ),
    default_cluster = "10x_graphclust",
    docs_path = "www"
)

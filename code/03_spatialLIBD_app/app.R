library("spatialLIBD")
library("markdown")
library("here")


# 2021-11-11T05:30:50.218127+00:00 shinyapps[5096402]: Warning: Error in loadNamespace: there is no package called ‘markdown’

## To install new spatialLIBD 1.15.4 (development version):
## .    https://bioconductor.org/packages/devel/data/experiment/html/spatialLIBD.html

## spatialLIBD uses golem.
## Golem is a framework for building production-grade shiny applications
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

here("code", "03_spatialLIBD_app")

## I added a symbolic link to point the spe.rds object to wrap.
spe <- readRDS("spe_harmony.rds") # spe with harmony and BayesSpace

# lobstr::obj_size(spe)
# 4.03 GB

## Import BayesSpace clusters
spe <- cluster_import(spe,
    cluster_dir = "clusters_BayesSpace",
    prefix = ""
)

## Quick inspection
# table(spe$BayesSpace_harmony_k24)

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
        "edge_spots",
        # vars[grep("^SNN_k10", vars)],
        # vars[grep("^BayesSpace_pca", vars)],
        vars[grep("^BayesSpace_harmony_", vars)]
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


## Note. If fails to read the rds object, go to Session Menu -> Set Working Directory -> To source File location

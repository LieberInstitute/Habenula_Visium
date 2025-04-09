library("spatialLIBD")
library("markdown")
library("here")

# CODE TO WRAP THE SPE WITH IN-TISSUE AND OUT-TISSUE SPOTS NONE FILTERED

## To install new spatialLIBD 1.15.4 (development version):
## .    https://bioconductor.org/packages/devel/data/experiment/html/spatialLIBD.html

## spatialLIBD uses golem.
## Golem is a framework for building production-grade shiny applications
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

here("code", "03_spatialLIBD_app_prefiltered")

## I added a symbolic link to point the spe*.rds object to wrap.
#spe <- readRDS(here('processed-data', '02_build_spe', "spe_raw.rds"))
spe <- readRDS("spe_scran_spotsweeper.rds")

## Sort samples to make them match with coding versions
lst_order <- sort(unique(spe$sample_id))
sample_order <- unlist(sapply(lst_order, function(i) {
    sort(unique(spe$sample_id)[grepl(i, unique(spe$sample_id))])
}))
sample_order

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
        vars[grep("^SNN_k10", vars)],
        vars[grep("*_outliers", vars)],
        # vars[grep("^BayesSpace_harmony_", vars)]
    ),
    spe_continuous_vars = c(
        "sum_umi",
        "sum_gene",
        "expr_chrM",
        "expr_chrM_ratio"
        # "edge_distance"
    ),
    default_cluster = "10x_graphclust",
    docs_path = "www"
)


## Note. If fails to read the rds object, go to Session Menu -> Set Working Directory -> To source File location

library("spatialLIBD")
library("markdown")
library("here")

## This the folder 03_spatialLIBD_app_deploy_k16/ to deploy the subset spatialLIBD shiny at BS k16

# 2021-11-11T05:30:50.218127+00:00 shinyapps[5096402]: Warning: Error in loadNamespace: there is no package called ‘markdown’

## To install new spatialLIBD 1.15.4 (development version):
## .    https://bioconductor.org/packages/devel/data/experiment/html/spatialLIBD.html

## spatialLIBD uses golem.
## Golem is a framework for building production-grade shiny applications
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Load the spe object
# spe_subset_for_spatialLIBD.rds
spe <- readRDS("spe_subset_for_spatialLIBD.rds")
# load the pseudobulked object sce_pseudo
sce_pseudo <- readRDS("sce_pseudo_BayesSpace_k16.rds")
# load modeling results for k9 clustering/pseudobulking
load("modeling_results_BayesSpace_k16.Rdata", verbose = TRUE)
load("sig_genes_subset_k16.Rdata", verbose = TRUE)
# lobstr::obj_size(spe)
# 4.03 GB

# Quick inspection
colnames(colData(spe))

## Import BayesSpace clusters
spe <- cluster_import(spe,
    cluster_dir = "clusters_BayesSpace",
    prefix = ""
)
spe$BayesSpace <- spe$BayesSpace_harmony_k16

## Quickly explore the data
vars <- colnames(colData(spe))
colnames(colData(spe)) <- vars <- gsub("X10x", "10x", vars)

colors_BayesSpace <- Polychrome::palette36.colors(28)
names(colors_BayesSpace) <- c(1:28)
m <- match(as.character(spe$BayesSpace_harmony_k16), names(colors_BayesSpace))
stopifnot(all(!is.na(m)))
spe$BayesSpace_colors <- spe$BayesSpace_harmony_k16_colors <- colors_BayesSpace[m]

spatialLIBD::run_app(
    spe,
    sce_layer = sce_pseudo,
    modeling_results = modeling_results,
    sig_genes = sig_genes,
    title = "spatialHabenula, Visium, Sp16",
    spe_discrete_vars = c( ## this is the variables for the spe object not the sce_pseudo object
        "BayesSpace",
        "ManualAnnotation",
        "overlaps_tissue",
        vars[grep("^10x_", vars)],
        vars[grep("^scran_", vars)],
        "edge_spots",
        vars[grep("^SNN_k10", vars)],
        # vars[grep("^BayesSpace_pca", vars)],
        vars[grep("^BayesSpace_harmony_", vars)],
        "BayesSpace_colors"
    ),
    spe_continuous_vars = c(
        "sum_umi",
        "sum_gene",
        "expr_chrM",
        "expr_chrM_ratio",
        "edge_distance"
    ),
    default_cluster = "BayesSpace", #"10x_graphclust",
    docs_path = "www"
)


## Note. If fails to read the rds object, go to Session Menu -> Set Working Directory -> To source File location

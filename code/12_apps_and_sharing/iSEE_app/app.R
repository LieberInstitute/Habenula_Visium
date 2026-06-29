library(SpatialExperiment)
library(iSEE)
library(shiny)
library(scuttle)
library(qs2)

#   For interactive testing at JHPCE
# setwd(here('code', '12_apps_and_sharing', 'iSEE_app'))

source("initial.R")

sce_pb = qs_read('sce_pb_shiny.qs2')

#   Use symbols for rownames but fall back on ENSEMBL for duplicates
rownames(sce_pb) <- uniquifyFeatureNames(
    rowData(sce_pb)$gene_id, rowData(sce_pb)$gene_name
)

## Don't run this on app.R since we don't want to run this every single time
# lobstr::obj_size(sce_pb)
# 20.83 MB

# sce_pb <- registerAppOptions(
#     sce_pb, color.maxlevels = length(colData(sce_pb)$BayesSpace_colors)
# )

iSEE(
    sce_pb,
    appTitle = "Habenula Atlas Pseudobulked Visium HD Data",
    initial = initial#,
    # colormap = ExperimentColorMap(colData = list(
    #     BayesSpace = function(n) {
    #         return(colData(sce_pb)$BayesSpace_colors)
    #     }
    # ))
)

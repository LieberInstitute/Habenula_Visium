library(SpatialExperiment)
library(iSEE)
library(shiny)
library(scuttle)
library(qs2)

#   For interactive testing at JHPCE
# setwd(here('code', '12_apps_and_sharing', 'iSEE_app'))

cell_type_colors = c(
    "MHb.1" = "#5e0c01",
    "MHb.2" = "#943f02",
    "Excit_LHb" = "#306171",
    "LHb.2.7" = "#ee9630",
    "LHb.4" = "#082844",
    "LHb.4/Inhib_LHb_4.2" = "#5e0c56",
    "Excit.Thal/Inhib_LHb_4.2" = "#8DADCA",
    "Excit.Thal" = "#2e6296",
    "Ependymal" = "#dbb369",
    "Subependymal" = "#ae9a7e",
    "Astrocyte" = "#972f2f",
    "Endo" = "#f65a45",
    "Endo/microglia" = "#141b02",
    "Oligo" =  "#384a08",
    "OPC" = "#829454"
)

source("initial.R")

sce_pb = qs_read('sce_pb_shiny.qs2')

#   Use symbols for rownames but fall back on ENSEMBL for duplicates
rownames(sce_pb) <- uniquifyFeatureNames(
    rowData(sce_pb)$gene_id, rowData(sce_pb)$gene_name
)

sce_pb$cell_type_colors = cell_type_colors[as.character(sce_pb$cell_type)]

## Don't run this on app.R since we don't want to run this every single time
# lobstr::obj_size(sce_pb)
# 20.83 MB

iSEE(
    sce_pb,
    appTitle = "Habenula Atlas Pseudobulked Visium HD Data",
    initial = initial,
    colormap = ExperimentColorMap(colData = list(
        cell_type = function(n) {
            return(sce_pb$cell_type_colors)
        }
    ))
)

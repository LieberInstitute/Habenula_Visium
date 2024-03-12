library("here")
library("spatialLIBD")
library("sessioninfo")


## output directory
dir_rdata <- here("processed-data", "04_harmony_BayesSpace")

## Load the data
spe <- readRDS(file.path(dir_rdata, "spe_harmony.rds"))

## Import BayesSpace clusters
spe <- cluster_import(spe,
    cluster_dir = file.path(dir_rdata, "clusters_BayesSpace"),
    prefix = ""
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

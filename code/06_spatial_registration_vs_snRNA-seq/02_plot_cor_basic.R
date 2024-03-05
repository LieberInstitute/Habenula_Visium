library("here")
library("spatialLIBD")
library("sessioninfo")

## Create output directories
dir_rdata <-
    here("processed-data", "06_spatial_registration_vs_snRNA-seq")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
dir_plots <-
    here("plots", "06_spatial_registration_vs_snRNA-seq")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

## Load correlation values
load(file.path(dir_rdata, "cor_BayesSpace_vs_snRNA-seq.Rdata"), verbose = TRUE)

## Make basic heatmaps (not ComplexHeatmap) versions

## Fine resolution
pdf(file = file.path(dir_plots, "snRNA-seq_registration_fineRes_basic.pdf"))
lapply(cor_fine, layer_stat_cor_plot, max = max(sapply(cor_fine, function(x) max(abs(x)))))
dev.off()

## Broad resolution
pdf(file = file.path(dir_plots, "snRNA-seq_registration_broadRes_basic.pdf"))
lapply(cor_broad, layer_stat_cor_plot, max = max(sapply(cor_broad, function(x) max(abs(x)))))
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

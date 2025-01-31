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
load(file.path(dir_rdata, "cor_BayesSpace_vs_snRNA-seq_pairwise.Rdata"),
    verbose = TRUE
)

## Make basic heatmaps (not ComplexHeatmap) versions

## Fine resolution
pdf(file = file.path(dir_plots, "snRNA-seq_registration_fineRes_basic_pairwise.pdf"))
lapply(
    cor_fine,
    layer_stat_cor_plot,
    max = max(sapply(cor_fine, max)),
    min = min(sapply(cor_fine, min))
)
dev.off()

## Broad resolution
pdf(file = file.path(dir_plots, "snRNA-seq_registration_broadRes_basic_pairwise.pdf"))
lapply(
    cor_broad,
    layer_stat_cor_plot,
    max = max(sapply(cor_broad, max)),
    min = min(sapply(cor_broad, min))
)
dev.off()


k08_broad <- cor_broad$BayesSpace_harmony_k08
rownames(k08_broad) <- gsub("^k08_", "", rownames(k08_broad))

pdf(file = file.path(dir_plots, "snRNA-seq_registration_broadRes_basic_k08_square_pairwise.pdf"))
layer_stat_cor_plot(k08_broad ,
    max = max(k08_broad),
    min = min(k08_broad))
dev.off()

pdf(file = file.path(dir_plots, "snRNA-seq_registration_broadRes_basic_k08_tall_pairwise.pdf"), height = 10)
layer_stat_cor_plot(k08_broad ,
    max = max(k08_broad),
    min = min(k08_broad))
dev.off()

pdf(file = file.path(dir_plots, "snRNA-seq_registration_broadRes_basic_k08_wide_pairwise.pdf"), width = 10)
layer_stat_cor_plot(k08_broad ,
    max = max(k08_broad),
    min = min(k08_broad))
dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

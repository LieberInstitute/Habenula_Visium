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

## Annotate clusters
annotated_clusters_fine <- lapply(cor_fine, annotate_registered_clusters, cutoff_merge_ratio = 0.1)
annotated_clusters_broad <- lapply(cor_broad, annotate_registered_clusters, cutoff_merge_ratio = 0.1)

## Use annotation labels on the correlation matrices
cor_fine_annotated <- mapply(function(cor, label_data) {
    rownames(cor) <-
        paste0(rownames(cor), " ~ ", label_data$layer_label[match(rownames(cor), label_data$cluster)])
    return(cor)
}, cor_fine, annotated_clusters_fine)

cor_broad_annotated <- mapply(function(cor, label_data) {
    rownames(cor) <-
        paste0(rownames(cor), " ~ ", label_data$layer_label[match(rownames(cor), label_data$cluster)])
    return(cor)
}, cor_broad, annotated_clusters_broad)



## Make basic heatmaps (not ComplexHeatmap) versions

## Fine resolution
pdf(file = file.path(dir_plots, "snRNA-seq_registration_fineRes_basic.pdf"))
lapply(cor_fine_annotated, layer_stat_cor_plot, max = max(sapply(cor_fine, max)), min = min(sapply(cor_fine, min)))
dev.off()

## Broad resolution
pdf(file = file.path(dir_plots, "snRNA-seq_registration_broadRes_basic.pdf"))
lapply(cor_broad_annotated, layer_stat_cor_plot, max = max(sapply(cor_broad, max)), min = min(sapply(cor_broad, min)))
dev.off()



## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

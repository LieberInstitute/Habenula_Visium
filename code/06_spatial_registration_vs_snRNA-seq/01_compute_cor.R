library("here")
library("spatialLIBD")
library("sessioninfo")

## Create output directories
dir_rdata <-
    here("processed-data", "06_spatial_registration_vs_snRNA-seq")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

## Load BayesSpace enrichment t-stats
BayesSpace_stats_list <-
    readRDS(here(
        "processed-data",
        "04_harmony_BayesSpace",
        "BayesSpace_stats_list.rds"
    ))

## Reformat to match expect format for spatialLIBD::layer_stat_cor
## Aka, match the format from
## head(spatialLIBD::tstats_Human_DLPFC_snRNAseq_Nguyen_topLayer)

BayesSpace_stats_list <- lapply(BayesSpace_stats_list, function(x) {
    rownames(x) <- x$gene
    x$gene <- NULL
    colnames(x) <- gsub("BayesSpace_harmony_", "", colnames(x))
    return(x)
})

registration_vars <-
    c("final_Annotations", "final_Annotations_broad")

compute_cor <- function(current_var) {
    ## Load input snRNA-seq data
    results_enrichment <-
        readRDS(here(
            "processed-data",
            "05_snRNA-seq_model_stats",
            paste0("enrichment_", current_var, ".rds")
        ))
    modeling_results_Hb <- list("enrichment" = results_enrichment)

    lapply(
        BayesSpace_stats_list,
        layer_stat_cor,
        modeling_results = modeling_results_Hb,
        top_n = 100
    )
}

cor_fine <- compute_cor("final_Annotations")
cor_broad <- compute_cor("final_Annotations_broad")

## Annotate clusters
annotated_clusters_fine <-
    lapply(cor_fine, annotate_registered_clusters, cutoff_merge_ratio = 0.1)
annotated_clusters_broad <-
    lapply(cor_broad, annotate_registered_clusters, cutoff_merge_ratio = 0.1)

## Use annotation labels on the correlation matrices
cor_fine <- mapply(function(cor, label_data) {
    rownames(cor) <-
        paste0(rownames(cor), " ~ ", label_data$layer_label[match(rownames(cor), label_data$cluster)])
    return(cor)
}, cor_fine, annotated_clusters_fine)

cor_broad <- mapply(function(cor, label_data) {
    rownames(cor) <-
        paste0(rownames(cor), " ~ ", label_data$layer_label[match(rownames(cor), label_data$cluster)])
    return(cor)
}, cor_broad, annotated_clusters_broad)

save(cor_fine,
    cor_broad,
    file = file.path(dir_rdata, "cor_BayesSpace_vs_snRNA-seq.Rdata")
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

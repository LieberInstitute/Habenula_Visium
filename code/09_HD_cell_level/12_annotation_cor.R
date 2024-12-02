library(here)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(HDF5Array)
library(sessioninfo)

t_stat_paths = here(
    'processed-data', '09_HD_cell_level', 'banksy', 'markers',
    sprintf('k%s.rds', 2:16)
)
out_path = here(
    'processed-data', '09_HD_cell_level', 'banksy', 'cor_vs_snRNA-seq.rds'
)
registration_vars = c("final_Annotations", "final_Annotations_broad")

compute_cor <- function(t_stats, current_var) {
    ## Load input snRNA-seq data
    results_enrichment <-
        readRDS(here(
            "processed-data",
            "05_snRNA-seq_model_stats",
            paste0("enrichment_", current_var, ".rds")
        ))
    modeling_results_Hb <- list("enrichment" = results_enrichment)

    lapply(
        t_stats,
        layer_stat_cor,
        modeling_results = modeling_results_Hb,
        top_n = 100
    )
}

t_stats = lapply(
    t_stat_paths,
    function(path) {
        readRDS(path) |>
            column_to_rownames('gene')
    }
)

cor_fine <- compute_cor(t_stats, registration_vars[1])
cor_broad <- compute_cor(t_stats, registration_vars[2])

## Annotate clusters
annotated_clusters_fine = lapply(
    cor_fine, annotate_registered_clusters, cutoff_merge_ratio = 0.1
)
annotated_clusters_broad = lapply(
    cor_broad, annotate_registered_clusters, cutoff_merge_ratio = 0.1
)

## Use annotation labels on the correlation matrices
cor_fine <- mapply(
    function(cor, label_data) {
        rownames(cor) = paste0(
            rownames(cor),
            " ~ ",
            label_data$layer_label[match(rownames(cor), label_data$cluster)]
        )
        return(cor)
    },
    cor_fine,
    annotated_clusters_fine
)

cor_broad <- mapply(
    function(cor, label_data) {
        rownames(cor) = paste0(
            rownames(cor),
            " ~ ",
            label_data$layer_label[match(rownames(cor), label_data$cluster)]
        )
        return(cor)
    },
    cor_broad,
    annotated_clusters_broad
)

cor_list = list(fine = cor_fine, broad = cor_broad)
saveRDS(cor_list, file = out_path)

session_info()

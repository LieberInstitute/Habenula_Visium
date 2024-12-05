library(here)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(HDF5Array)
library(sessioninfo)

t_stat_paths = here(
    'processed-data', '09_HD_cell_level', 'banksy', 'markers',
    sprintf('k%s.rds', 2:28)
)
out_path = here(
    'processed-data', '09_HD_cell_level', 'banksy', 'cor_vs_snRNA-seq.rds'
)
plot_dir = here('plots', '09_HD_cell_level', 'banksy')
registration_vars = c("final_Annotations", "final_Annotations_broad")
res_names = c('fine', 'broad')

annotated_heatmap <- function(t_stats, current_var, res_name) {
    ## Load input snRNA-seq data
    results_enrichment <-
        readRDS(here(
            "processed-data",
            "05_snRNA-seq_model_stats",
            paste0("enrichment_", current_var, ".rds")
        ))
    modeling_results_Hb <- list("enrichment" = results_enrichment)

    this_cor = lapply(
        t_stats,
        layer_stat_cor,
        modeling_results = modeling_results_Hb,
        top_n = 100
    )

    #   Put clusters in order
    this_cor = lapply(
        this_cor,
        function(x) {
            x[
                rownames(x) |>
                    str_extract('_([0-9]+)$', group = 1) |>
                    as.numeric() |>
                    order(),
            ]
        }
    )

    #   Annotate clusters
    annotated_clusters = lapply(
        this_cor, annotate_registered_clusters, cutoff_merge_ratio = 0.1
    )

    #   Use annotation labels on the correlation matrices
    this_cor <- mapply(
        function(cor, label_data) {
            rownames(cor) = paste0(
                rownames(cor),
                " ~ ",
                label_data$layer_label[match(rownames(cor), label_data$cluster)]
            )
            return(cor)
        },
        this_cor,
        annotated_clusters
    )

    #   Make basic heatmap (not ComplexHeatmap)
    pdf(
        file.path(
            plot_dir,
            sprintf("snRNA-seq_registration_%sRes_basic.pdf", res_name)
        )
    )
    lapply(
        this_cor,
        layer_stat_cor_plot,
        max = max(sapply(this_cor, max)),
        min = min(sapply(this_cor, min))
    )
    dev.off()

    return(this_cor)
}

#   Read in t stats for genes marking Banksy clusters at each value of k
t_stats = lapply(
    t_stat_paths,
    function(path) {
        readRDS(path) |>
            column_to_rownames('gene')
    }
)

#   Compute and annotate correlation matrices comparing t stats of Banksy
#   clusters against the snRNA-seq data. Save heatmaps and the matrices
#   themselves
cor_fine = annotated_heatmap(t_stats, registration_vars[1], res_names[1])
cor_broad = annotated_heatmap(t_stats, registration_vars[2], res_names[2])
cor_list = list(fine = cor_fine, broad = cor_broad)
saveRDS(cor_list, file = out_path)

session_info()

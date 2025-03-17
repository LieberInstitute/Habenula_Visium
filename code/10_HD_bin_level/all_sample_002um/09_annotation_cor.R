

#####
#remotes::install_github('LieberInstitute/spatialLIBD')

library(here)
library(tidyverse)
library(SpatialExperiment)
library(HDF5Array)
library(sessioninfo)
library(spatialLIBD)

t_stat_paths = here('processed-data', '09_HD_cell_level', 'ficture', 'markers', 'all_samples_12.rds')
out_path = here('processed-data', '09_HD_cell_level', 'ficture_aggregate', "cor_vs_snRNA-seq_all_samples.rds")

plot_dir = here('plots', '09_HD_cell_level', 'ficture_aggregate')
plot_dir2 = here('plots', '09_HD_cell_level', 'ficture_aggregate','1')
registration_vars = c("final_Annotations", "final_Annotations_broad")
res_names = c('fine', 'broad')
dir.create(dirname(out_path), showWarnings = FALSE)
dir.create(dirname(plot_dir2), showWarnings = FALSE)

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
    # this_cor <- mapply(
    # function(cor, label_data) {
    #     rownames(cor) <- paste0(
    #         rownames(cor),
    #         " ~ ",
    #         label_data$layer_label[match(rownames(cor), label_data$cluster)]
    #     )
    #     return(cor)
    # },
    # list(this_cor[[1]]),
    # list(annotated_clusters[[1]]),
    # SIMPLIFY = FALSE
    # )
    
    #   Make basic heatmap (not ComplexHeatmap)
    pdf(
        file.path(
            plot_dir,
            sprintf("snRNA-seq_registration_%sRes_basic.pdf", res_name)
        )
    )

    print(layer_stat_cor_plot(this_cor[[1]], annotation = annotated_clusters[[1]],cluster_rows = FALSE,
    cluster_columns = FALSE))
    

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

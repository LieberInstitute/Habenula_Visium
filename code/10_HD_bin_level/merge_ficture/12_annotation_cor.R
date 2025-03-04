library(here)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(HDF5Array)
library(sessioninfo)

task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
sample_id<-c("H1-W369TJK_D1_9090","H1-MVPY9BW_A1_8433","H1-MVPY9BW_D1_8667","H1-XQQD7C7_A1_8518","H1-XQQD7C7_D1_9037")
i<- sample_id[task_id]

t_stat_paths = here('processed-data', '10_HD_bin_level', 'ficture', 'markers', sprintf('%s_12.rds', as.character(i)))
out_path = here('processed-data', '10_HD_bin_level', 'banksy', sprintf('cor_vs_snRNA-seq_%s.rds', as.character(i)))
plot_dir = here('plots', '10_HD_bin_level', 'banksy')
plot_dir2 = here('plots', '10_HD_bin_level', 'banksy','1')
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
    this_cor <- mapply(
    function(cor, label_data) {
        rownames(cor) <- paste0(
            rownames(cor),
            " ~ ",
            label_data$layer_label[match(rownames(cor), label_data$cluster)]
        )
        return(cor)
    },
    list(this_cor[[1]]),
    list(annotated_clusters[[1]]),
    SIMPLIFY = FALSE
    )

    #   Make basic heatmap (not ComplexHeatmap)
    pdf(
        file.path(
            plot_dir,
            sprintf("snRNA-seq_registration_%sRes_basic_%s.pdf", res_name,i)
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
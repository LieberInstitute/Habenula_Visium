library(here)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(sessioninfo)

model_paths = here(
    'processed-data', '09_HD_cell_level', 'registration_banksy',
    'modeling_results',
    sprintf('%s.rds', sub('\\.', '_', as.character(seq_len(10) / 10)))
)
ref_paths = here(
    "processed-data", "05_snRNA-seq_model_stats",
    sprintf(
        "enrichment_%s.rds", c("final_Annotations", "final_Annotations_broad")
    )
)
out_path = here(
    'processed-data', '09_HD_cell_level', 'registration_banksy',
    'cor_vs_snRNA-seq.rds'
)
plot_dir = here('plots', '09_HD_cell_level', 'registration_banksy')
res_names = c('fine', 'broad')

dir.create(plot_dir, showWarnings = FALSE)

annotated_heatmap <- function(t_stats, ref_path, res_name) {
    ## Load input snRNA-seq data
    results_enrichment <- readRDS(ref_path)
    modeling_results_sn <- list("enrichment" = results_enrichment)

    this_cor = lapply(
        t_stats,
        layer_stat_cor,
        modeling_results = modeling_results_sn,
        model_type = "enrichment",
        top_n = 100
    )

    #   Annotate clusters
    annotated_clusters = lapply(
        this_cor, annotate_registered_clusters, cutoff_merge_ratio = 0.1
    )

    #   Make heatmaps
    pdf(
        file.path(
            plot_dir,
            sprintf("snRNA-seq_registration_%sRes_basic.pdf", res_name)
        )
    )
    for (i in seq_len(length(this_cor))) {
        print(layer_stat_cor_plot(this_cor[[i]], annotation = annotated_clusters[[i]]))
    }
    dev.off()

    return(this_cor)
}

#   Read in enrichment stats for Banksy clusters at all Leiden resolutions
t_stats = lapply(model_paths, function(path) readRDS(path)$enrichment)

#   Compute and annotate correlation matrices comparing t stats of Banksy
#   clusters against the snRNA-seq data. Save heatmaps and the matrices
#   themselves
cor_fine = annotated_heatmap(t_stats, ref_paths[1], res_names[1])
cor_broad = annotated_heatmap(t_stats, ref_paths[2], res_names[2])
cor_list = list(fine = cor_fine, broad = cor_broad)
saveRDS(cor_list, file = out_path)

session_info()

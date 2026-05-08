#   Correlation heatmaps: registration of each intracellular vs extracellular
#   FICTURE result

library(here)
library(tidyverse)
library(spatialLIBD)
library(sessioninfo)

k = c(3:10, 20)

plot_dir = here(
    'plots', '10_HD_bin_level', 'no_secondary', 'cell_environment', 
    'investigation', 'registration'
)
model_paths = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'investigation', 'registration', sprintf('%s_modeling.rds', k)
)
ref_paths = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'registration', 'modeling_results', sprintf('%s.rds', k)
)
out_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'investigation', 'registration', 'cor_results.rds'
)

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(out_path), showWarnings = FALSE)

#   Read in enrichment stats for FICTURE clusters at all k values, as well as
#   the reference data
t_stats = lapply(model_paths, function(path) readRDS(path)$enrichment)
results_enrichment = lapply(
    ref_paths, function(path) list(enrichment = readRDS(path)$enrichment)
)

this_cor = map2(
    t_stats, results_enrichment,
    ~ layer_stat_cor(.x, modeling_results = .y, model_type = "enrichment")
)

#  Clean up cluster names
for (i in seq_len(length(this_cor))) {
    rownames(this_cor[[i]]) = sub('^X', 'Factor_', rownames(this_cor[[i]]))
    colnames(this_cor[[i]]) = sub('^X', 'Factor_', colnames(this_cor[[i]]))

}

#   Annotate clusters
annotated_clusters = lapply(
    this_cor, annotate_registered_clusters, cutoff_merge_ratio = 0.1
)

#   Make heatmaps
pdf(file.path(plot_dir, "cor_heatmaps.pdf"))
for (i in seq_len(length(this_cor))) {
    print(
        layer_stat_cor_plot(
            this_cor[[i]], annotation = annotated_clusters[[i]]#,
            # heatmap_legend_param = list(title = "Cor", at = c(-1, 0, 1))
        )
    )
}
dev.off()

saveRDS(this_cor, file = out_path)

session_info()

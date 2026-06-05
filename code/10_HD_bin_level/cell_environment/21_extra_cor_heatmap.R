#   Correlation heatmaps: registration for each FICTURE k value against Banksy
#   cell types (not clusters)

library(here)
library(tidyverse)
library(spatialLIBD)
library(sessioninfo)

k = 3:20

plot_dir = here(
    'plots', '10_HD_bin_level', 'no_secondary', 'cell_environment', 
    'registration'
)
model_paths = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'registration', 'modeling_results', sprintf('%s.rds', k)
)
ref_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'modeling_results', '1_8_cell_types.rds'
)
out_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'registration', 'cor_vs_cell_types.rds'
)

dir.create(plot_dir, showWarnings = FALSE)
dir.create(dirname(out_path), showWarnings = FALSE)

#   Read in enrichment stats for FICTURE clusters at all k values, as well as
#   the reference data
t_stats = lapply(model_paths, function(path) readRDS(path)$enrichment)
results_enrichment = list(enrichment = readRDS(ref_path)$enrichment)

this_cor = lapply(
    t_stats,
    layer_stat_cor,
    modeling_results = results_enrichment,
    model_type = "enrichment"
)

#  Remove 'X' from cluster names
for (i in seq_len(length(this_cor))) {
    rownames(this_cor[[i]]) = sub('^X', '', rownames(this_cor[[i]]))
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
            this_cor[[i]], annotation = annotated_clusters[[i]],
            heatmap_legend_param = list(title = "Cor", at = c(-1, 0, 1))
        )
    )
}
dev.off()

saveRDS(this_cor, file = out_path)

session_info()

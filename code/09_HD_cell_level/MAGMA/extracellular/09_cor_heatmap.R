#   Generate correlation heatmap for extracellular vs cellular cell types

library(here)
library(tidyverse)
library(spatialLIBD)
library(sessioninfo)

plot_dir = here(
    'plots', '09_HD_cell_level', 'no_secondary', 'MAGMA',
    'extracellular'
)
model_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA',
    'extracellular', 'spatial_registration', 'modeling_results.rds'
)
ref_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'modeling_results', '1_8_cell_types.rds'
)
out_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA',
    'extracellular', 'spatial_registration', 'cor_vs_cellular.rds'
)

dir.create(plot_dir, showWarnings = FALSE)
dir.create(dirname(out_path), showWarnings = FALSE)

t_stats_target = readRDS(model_path)$enrichment
t_stats_ref = list(enrichment = readRDS(ref_path)$enrichment)

this_cor = layer_stat_cor(
    t_stats_target, modeling_results = t_stats_ref, model_type = "enrichment"
)

anno_df = annotate_registered_clusters(this_cor, cutoff_merge_ratio = 0.1)

#   Make heatmap
pdf(file.path(plot_dir, "cor_heatmap.pdf"))
print(
    layer_stat_cor_plot(
        this_cor, annotation = anno_df,
        heatmap_legend_param = list(title = "Cor", at = c(-1, 0, 1)),
        row_title = "Extracellular Cell Type",
        column_title = "Cellular Cell Type"
    )
)
dev.off()

saveRDS(this_cor, file = out_path)

session_info()

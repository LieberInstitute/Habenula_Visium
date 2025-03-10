library(here)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(sessioninfo)

plot_dir = here('plots', '09_HD_cell_level', 'registration_banksy')
model_paths = here(
    'processed-data', '09_HD_cell_level', 'registration_banksy',
    'modeling_results',
    sprintf('%s.rds', sub('\\.', '_', as.character(seq_len(10) / 10)))
)

#   List all paths and names for reference data
ref_paths = here(
    "processed-data", "05_snRNA-seq_model_stats",
    sprintf(
        "enrichment_%s.rds", c("final_Annotations", "final_Annotations_broad")
    )
)
ref_names = c('snRNAseq_fine', 'snRNAseq_broad')

#   Get the reference data for this task
task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
ref_path = ref_paths[task_id]
ref_name = ref_names[task_id]

out_path = here(
    'processed-data', '09_HD_cell_level', 'registration_banksy',
    sprintf('cor_vs_%s.rds', ref_name)
)

dir.create(plot_dir, showWarnings = FALSE)

#   Read in enrichment stats for Banksy clusters at all Leiden resolutions
t_stats = lapply(model_paths, function(path) readRDS(path)$enrichment)

#   Load reference data
results_enrichment <- readRDS(ref_path)
modeling_results_sn <- list("enrichment" = results_enrichment)

#   Correlate Banksy clusters with reference data
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
pdf(file.path(plot_dir, sprintf("%s.pdf", ref_name)))
for (i in seq_len(length(this_cor))) {
    print(layer_stat_cor_plot(this_cor[[i]], annotation = annotated_clusters[[i]]))
}
dev.off()

saveRDS(this_cor, file = out_path)

session_info()

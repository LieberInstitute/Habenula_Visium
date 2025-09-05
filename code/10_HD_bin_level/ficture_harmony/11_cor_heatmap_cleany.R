library(here)
library(tidyverse)
library(SpatialExperiment)
library(HDF5Array)
library(sessioninfo)
library(spatialLIBD)

model_paths = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'ficture_harmony',
    'registration', 'modeling_results', 'cleaning_y',
    sprintf('%d.rds', c(seq(3, 40), 70, 100))
)
plot_dir = here(
    'plots', '10_HD_bin_level', 'probe_fix', 'ficture_harmony', 'registration'
)

#   List all paths and names for reference data
ref_paths = c(
    #   Fine and broad snRNA-seq data
    here(
        "processed-data", "05_snRNA-seq_model_stats",
        sprintf(
            "enrichment_%s.rds",
            c("final_Annotations", "final_Annotations_broad")
        )
    ),
    #   Multiome data
    here(
        'processed-data', '05_snRNA-seq_model_stats',
        'enrichment_snRNA-multiome_v2.rds'
    ),
    #    Visium BayesSpace clusters (k 2 through 28)
    here(
        "processed-data", "05_brain_area_differential_expression",
        "modeling_results_BS",
        sprintf("modeling_results_BayesSpace_k%02d.Rdata", 2:28)
    )
)
ref_names = c(
    'snRNAseq_fine', 'snRNAseq_broad', 'multiome',
    sprintf('Visium_BayesSpace_k%02d', 2:28)
)

#   Get the reference data for this task
task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
ref_path = ref_paths[task_id]
ref_name = ref_names[task_id]

out_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'ficture_harmony',
    'registration', 'cor_rds', 'cleaning_y', sprintf('cor_vs_%s.rds', ref_name)
)

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

#   Read in enrichment stats for Banksy clusters at all Leiden resolutions
t_stats = lapply(model_paths, function(path) readRDS(path)$enrichment)

#   Load reference data
if (grepl('^Visium', ref_name)) {
    results_enrichment = get(load(ref_path))
} else if (grepl('^snRNAseq', ref_name)) {
    results_enrichment = list(enrichment = readRDS(ref_path))
} else {
    results_enrichment = list(
        enrichment = readRDS(ref_path)$enrichment |>
            filter(!duplicated(ensembl))
    )
}

# message("Checking t_stats[[1]] structure...")
# print(str(t_stats[[1]]))

# message("Checking results_enrichment$enrichment structure...")
# print(str(results_enrichment$enrichment))

# message("R options:")
# print(options("stringsAsFactors"))

# message("R locale:")
# print(Sys.getlocale())

# message("Package version:")
# print(packageVersion("spatialLIBD"))

this_cor = lapply(
    t_stats,
    layer_stat_cor,
    modeling_results = results_enrichment,
    model_type = "enrichment",
    top_n = 100
)

#  Remove 'X' from cluster names
for (i in seq_len(length(this_cor))) {
    rownames(this_cor[[i]]) = sub('^X', '', rownames(this_cor[[i]]))
}

for (i in seq_len(length(this_cor))) {
    rownames(this_cor[[i]]) = sub('^c', '', rownames(this_cor[[i]]))
}

#   Annotate clusters
annotated_clusters = lapply(
    this_cor, annotate_registered_clusters, cutoff_merge_ratio = 0.1
)

#   Make heatmaps
pdf(file.path(plot_dir, sprintf("%s.pdf", ref_name)))
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

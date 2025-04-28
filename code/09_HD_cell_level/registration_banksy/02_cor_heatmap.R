library(here)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(sessioninfo)
library(getopt)

# Import command-line parameters
spec <- matrix(
    c(
        c("ref", "lambda"),
        c("r", "l"),
        rep("1", 2),
        c("integer", "numeric"),
        c("Reference data index", "Banksy lambda parameter")
    ),
    ncol = 5
)
opt <- getopt(spec)

message("Using the following parameters:")
print(opt)

lambda_neat = paste0('lambda', sub('\\.', '_', as.character(opt$lambda)))

plot_dir = here(
    'plots', '09_HD_cell_level', 'probe_fix', 'registration_banksy', lambda_neat
)
model_paths = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'registration_banksy',
    'modeling_results', lambda_neat,
    sprintf('%s.rds', c(sub('\\.', '_', as.character(seq_len(20) / 10)), 4, 8))
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
        "processed-data", "05_layer_differential_expression",
        "modeling_results_BS",
        sprintf("modeling_results_BayesSpace_k%02d.Rdata", 2:28)
    )
)
ref_names = c(
    'snRNAseq_fine', 'snRNAseq_broad', 'multiome',
    sprintf('Visium_BayesSpace_k%02d', 2:28)
)

#   Get the reference data for this task
ref_path = ref_paths[opt$ref]
ref_name = ref_names[opt$ref]

out_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'registration_banksy',
    lambda_neat, sprintf('cor_vs_%s.rds', ref_name)
)

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(out_path), showWarnings = FALSE)

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

#   For a talk for Kristen, produce custom correlation heatmaps for a specific
#   Banksy resolution, just the multiome and fine snRNA-seq references, and
#   only habenula- and thalamus-related clusters

library(here)
library(tidyverse)
library(spatialLIBD)
library(sessioninfo)

plot_dir = here(
    'plots', '09_HD_cell_level', 'probe_fix', 'registration_banksy', 'lambda0_2'
)
model_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'registration_banksy',
    'modeling_results', 'lambda0_2', '1_4.rds'
)

#   List all paths and names for reference data
ref_paths = c(
    #   Fine snRNA-seq data
    here(
        "processed-data", "05_snRNA-seq_model_stats",
        "enrichment_final_Annotations.rds"
    ),
    #   Multiome data
    here(
        'processed-data', '05_snRNA-seq_model_stats',
        'enrichment_snRNA-multiome_v2.rds'
    )
)
ref_names = c('snRNAseq_fine', 'multiome')

#   Get the reference data for this task
index = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
ref_path = ref_paths[index]
ref_name = ref_names[index]

#   Read in enrichment stats for Banksy clusters
t_stats = readRDS(model_path)$enrichment

#   Load reference data
if (ref_name == 'snRNAseq_fine') {
    results_enrichment = list(enrichment = readRDS(ref_path))
} else {
    results_enrichment = list(
        enrichment = readRDS(ref_path)$enrichment |>
            filter(!duplicated(ensembl))
    )
}

this_cor = layer_stat_cor(
    t_stats,
    modeling_results = results_enrichment,
    model_type = "enrichment",
    top_n = 100
)

#  Remove 'X' from cluster names
rownames(this_cor) = sub('^X', '', rownames(this_cor))

#   Annotate clusters
annotated_clusters = annotate_registered_clusters(
    this_cor, cutoff_merge_ratio = 0.1
)

if (ref_name == 'snRNAseq_fine') {
    #   Filter to habenula and thalamus clusters only, then take only clusters
    #   with at least one X
    annotated_clusters = annotated_clusters |>
        filter(
            grepl("[ML]Hb|\\.Thal", layer_label),
            layer_confidence == 'good'
        )
    this_cor = this_cor |>
        as.data.frame() |>
        rownames_to_column('cluster') |>
        select(cluster, matches("^[ML]Hb|\\.Thal$")) |>
        filter(cluster %in% annotated_clusters$cluster) |>
        column_to_rownames('cluster') |>
        as.matrix()
} else {
    #   Filter to habenula clusters only, then take only clusters with at least
    #   one X
    annotated_clusters = annotated_clusters |>
        filter(
            grepl("DD_[ML]Hb", layer_label),
            layer_confidence == 'good'
        )
    this_cor = this_cor |>
        as.data.frame() |>
        rownames_to_column('cluster') |>
        select(cluster, matches("DD_[ML]Hb$")) |>
        filter(cluster %in% annotated_clusters$cluster) |>
        column_to_rownames('cluster') |>
        as.matrix()
}

#   Make heatmaps
pdf(file.path(plot_dir, sprintf("custom_%s.pdf", ref_name)))
layer_stat_cor_plot(
    this_cor, annotation = annotated_clusters,
    heatmap_legend_param = list(title = "Cor", at = c(-1, 0, 1))
)
dev.off()

session_info()

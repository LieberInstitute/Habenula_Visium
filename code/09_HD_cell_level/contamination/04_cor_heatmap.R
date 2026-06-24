#   Generate correlation heatmaps for spatial registration against snRNA-seq,
#   multiome, and Visium BayesSpace reference datasets

library(here)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(sessioninfo)

plot_dir = here(
    'plots', '09_HD_cell_level', 'no_secondary', 'contamination'
)
model_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'contamination',
    'registration', 'modeling_results.rds'
)
ref_path = '/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/11_link_prep/04_registration_wrapper/model_results_mid.rds'
out_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'contamination',
    'registration', 'cor.rds'
)

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(out_path), showWarnings = FALSE)

t_stats = readRDS(model_path)$enrichment

#   This duplicate-gene-id thing is a bug and should be fixed, but it only
#   affects 9 genes, is complex to fix, and not a priority for this analysis
results_enrichment = list(
    enrichment = readRDS(ref_path)$enrichment |>
        filter(!duplicated(ensembl))
)

this_cor = layer_stat_cor(
    t_stats, modeling_results = results_enrichment, model_type = "enrichment"
)
rownames(this_cor) = sub('^X', '', rownames(this_cor))

annotated_clusters = annotate_registered_clusters(
    this_cor, cutoff_merge_ratio = 0.1
)

pdf(file.path(plot_dir, "full_heatmap.pdf"))
print(
    layer_stat_cor_plot(
        this_cor, annotation = annotated_clusters,
        heatmap_legend_param = list(title = "Cor", at = c(-1, 0, 1))
    )
)
dev.off()

pdf(file.path(plot_dir, "11_only_heatmap.pdf"), height = 4)
print(
    layer_stat_cor_plot(
        this_cor[grepl('^11_', rownames(this_cor)),],
        annotation = annotated_clusters[
            grepl('^11_', annotated_clusters$cluster),
        ]
    )
)
dev.off()


saveRDS(this_cor, file = out_path)

session_info()

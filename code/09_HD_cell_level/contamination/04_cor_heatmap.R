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
ct_anno_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'cluster_annotation.csv'
)
ref_path = '/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/11_link_prep/04_registration_wrapper/model_results_mid.rds'
out_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'contamination',
    'registration', 'cor.rds'
)
cell_type_colors = c(
    'Excit.Thal' = "#4d55b7",
    'Excit.Thal/Inhib_LHb_4.2' = "#2d1e6a",
    'LHb.4' = "#00607A",
    'LHb.4/Inhib_LHb_4.2' = "#003d4e",
    'Inhib.Thal' = "#9a9fe7",
    'Astrocyte' = "#532222",
    'OPC' = "#829454",
    'Oligo' =  "#384a08",
    'Microglia' = "#141b02",
    'LHb.2.7' = "#6C9FA9",
    'Endo' = "#d95f02",
    'Endo/microglia' = "#8d3e01",
    'Excit_LHb' = "#A8B8BC",
    'MHb.1' = "#A86A9A",
    'MHb.2' = "#BCA6B6",
    'Ependymal' = "#f5a105ff",
    'Subependymal' = "#976d1d"
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

#   Normal spatial registration heatmaps-- all clusters then just focusing on
#   cluster 11
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

#   As a crude indirect measure of possible contamination, compare nuclear and
#   cellular t-stats for each cluster via correlation
p = t_stats |>
    as_tibble() |>
    select(ensembl, matches("^t_stat_")) |>
    pivot_longer(
        cols = matches("^t_stat_"),
        names_to = c("cluster", "compartment"),
        names_pattern = "t_stat_X(\\d+)_(\\w+)",
        values_to = "t_stat"
    ) |>
    mutate(cluster = as.integer(cluster)) |>
    pivot_wider(
        names_from = compartment, values_from = t_stat, names_prefix = "t_stat_"
    ) |>
    group_by(cluster) |>
    summarize(
        cor_val = cor(t_stat_cellular, t_stat_nuclear, method = "pearson")
    ) |>
    left_join(read_csv(ct_anno_path, show_col_types = FALSE), by = 'cluster') |>
    arrange(desc(cor_val)) |>
    mutate(
        cluster = sprintf('%s ~ %s', fine_cell_type, cluster),
        cluster = factor(cluster, levels = cluster)
    ) |>
    ggplot(aes(x = cluster, y = cor_val, fill = fine_cell_type)) +
        geom_col() +
        scale_fill_manual(values = cell_type_colors) +
        guides(fill = "none") +
        theme_bw(base_size = 15) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
pdf(
    file.path(plot_dir, "cellular_nuclear_correlation.pdf"),
    height = 6, width = 8
)
print(p)
dev.off()

saveRDS(this_cor, file = out_path)

session_info()

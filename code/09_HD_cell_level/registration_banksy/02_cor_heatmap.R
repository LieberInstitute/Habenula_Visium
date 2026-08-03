#   Generate correlation heatmaps for spatial registration against snRNA-seq,
#   multiome, and Visium BayesSpace reference datasets

library(here)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(sessioninfo)

all_genes = TRUE

if (all_genes) {
    plot_dir = here(
        'plots', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
        'all_genes'
    )
} else {
    plot_dir = here(
        'plots', '09_HD_cell_level', 'no_secondary', 'registration_banksy'
    )
}

model_paths = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'modeling_results',
    sprintf(
        '%s.rds',
        c(sub('\\.', '_', as.character(seq_len(20) / 10)), 4, 8)
    )
)
cell_map_path = '/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/raw-data/cell_type_map.csv'
hd_cell_map_path = here('raw-data', 'sample_info', 'hd_cell_type_map.csv')
anno_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'cluster_annotation.csv'
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
    #   Multiome data mid and fine resolutions
    sprintf(
        '/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/11_link_prep/04_registration_wrapper/model_results_%s.rds',
        c("mid", "fine")
    ),
    #    Visium BayesSpace clusters (k 2 through 28)
    here(
        "processed-data", "05_brain_area_differential_expression",
        "modeling_results_BS",
        sprintf("modeling_results_BayesSpace_k%02d.Rdata", 2:28)
    )
)
ref_names = c(
    'snRNAseq_fine', 'snRNAseq_broad', 'multiome_mid', 'multiome_fine',
    sprintf('Visium_BayesSpace_k%02d', 2:28)
)

#   Get the reference data for this task
array_task = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
ref_path = ref_paths[array_task]
ref_name = ref_names[array_task]

if (all_genes) {
    out_path = here(
        'processed-data', '09_HD_cell_level', 'no_secondary',
        'registration_banksy', 'all_genes', sprintf('cor_vs_%s.rds', ref_name)
    )
} else {
    out_path = here(
        'processed-data', '09_HD_cell_level', 'no_secondary',
        'registration_banksy', sprintf('cor_vs_%s.rds', ref_name)
    )
}

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

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

if (all_genes) {
    this_cor = lapply(
        t_stats,
        layer_stat_cor,
        modeling_results = results_enrichment,
        model_type = "enrichment"
    )
} else {
    this_cor = lapply(
        t_stats,
        layer_stat_cor,
        modeling_results = results_enrichment,
        model_type = "enrichment",
        top_n = 100
    )
}

#  Remove 'X' from cluster names
for (i in seq_len(length(this_cor))) {
    rownames(this_cor[[i]]) = sub('^X', '', rownames(this_cor[[i]]))
}

#   Cell-type names were updated for multiome mid reference. Also, for the
#   manuscript, use better names for HD clusters
if (ref_name == 'multiome_mid') {
    cell_map_df = read_csv(cell_map_path, show_col_types = FALSE)
    anno_df = read_csv(anno_path, show_col_types = FALSE) |>
        left_join(
            read_csv(hd_cell_map_path, show_col_types = FALSE),
            by = c('fine_cell_type' = 'old_cell_type')
        )

    this_cor_fancy = lapply(
        this_cor,
        function(x) {
            colnames(x) = tibble(old_cell_type = colnames(x)) |>
                left_join(cell_map_df, by = 'old_cell_type') |>
                pull(new_cell_type)
            rownames(x) = tibble(cluster = as.numeric(rownames(x))) |>
                left_join(anno_df, by = 'cluster') |>
                mutate(
                    cluster_label = sprintf(
                        '%02d ~ %s', cluster, new_cell_type
                    )
                ) |>
                pull(cluster_label)
            return(x)
        }
    )

    #   Annotate clusters
    annotated_clusters = lapply(
        this_cor_fancy, annotate_registered_clusters, cutoff_merge_ratio = 0.1
    )

    #   Make heatmaps
    pdf(file.path(plot_dir, sprintf("%s_manuscript.pdf", ref_name)))
    for (i in seq_len(length(this_cor_fancy))) {
        print(
            layer_stat_cor_plot(
                this_cor_fancy[[i]], annotation = annotated_clusters[[i]]
            )
        )
    }
    dev.off()
} else if (ref_name == 'multiome_fine') {
    anno_df = read_csv(anno_path, show_col_types = FALSE) |>
        left_join(
            read_csv(hd_cell_map_path, show_col_types = FALSE),
            by = c('fine_cell_type' = 'old_cell_type')
        )

    this_cor_fancy = lapply(
        this_cor,
        function(x) {
            rownames(x) = tibble(cluster = as.numeric(rownames(x))) |>
                left_join(anno_df, by = 'cluster') |>
                mutate(
                    cluster_label = sprintf(
                        '%02d ~ %s', cluster, new_cell_type
                    )
                ) |>
                pull(cluster_label)
            return(x)
        }
    )

    #   Annotate clusters
    annotated_clusters = lapply(
        this_cor_fancy, annotate_registered_clusters, cutoff_merge_ratio = 0.1
    )

    #   Make heatmaps
    pdf(file.path(plot_dir, sprintf("%s_manuscript.pdf", ref_name)))
    for (i in seq_len(length(this_cor_fancy))) {
        print(
            layer_stat_cor_plot(
                this_cor_fancy[[i]], annotation = annotated_clusters[[i]]
            )
        )
    }
    dev.off()
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

#   For helping to annotate the multiome data, we also want a version where
#   Visium HD is the reference
if (grepl('^multiome', ref_name)) {
    #   Remove 'X' from Visium HD cluster names
    for (i in seq_len(length(t_stats))) {
        colnames(t_stats[[i]]) = sub(
            '_X([0-9]+)', '_\\1', colnames(t_stats[[i]])
        )
    }

    this_cor = lapply(
        t_stats,
        function(x) {
            layer_stat_cor(
                results_enrichment$enrichment,
                modeling_results = list(enrichment = x),
                model_type = "enrichment",
                top_n = 100
            )
        }
    )

    #   Annotate clusters
    annotated_clusters = lapply(
        this_cor, annotate_registered_clusters, cutoff_merge_ratio = 0.1
    )

    #   Make heatmaps
    pdf(file.path(plot_dir, sprintf("%s_flipped.pdf", ref_name)))
    for (i in seq_len(length(this_cor))) {
        print(
            layer_stat_cor_plot(
                this_cor[[i]], annotation = annotated_clusters[[i]],
                heatmap_legend_param = list(title = "Cor", at = c(-1, 0, 1))
            )
        )
    }
    dev.off()
}

session_info()

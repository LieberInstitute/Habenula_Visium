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

visium_manual_row_order = sprintf('Sp09D0%s', c(4, 8, 6, 2, 9, 5, 3, 7, 1))
visium_manual_col_order = as.character(
    c(6, 1, 14, 16, 11, 4, 5, 9, 13, 15, 3, 12, 17, 18, 10, 7, 2, 8)
)
multiome_manual_row_order = paste0(
    'C.',
    c(
        '16.DD_MHb', '11.DD_MHb', '07.DD_MHb', '10.DD_MHb', '14.DD_MHb', '34',
        '02', '22', '29', '26', '20', '21', '27', '41', '05.DD_LHb',
        '18.DD_LHb', '23.DD_LHb', '33.DD_LHb', '08', '36.DD_MHb', '24.DD_LHb',
        '30.DD_LHb', '40.DD_LHb', '13', '09', '31', '06', '04', '32', '12',
        '25', '17', '15', '35', '03', '37', '39', '28', '38', '19', '01'
    )
)

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
t_stats = lapply(model_paths, function(path) readRDS(path))

#   Load reference data
if (grepl('^Visium', ref_name)) {
    results_enrichment = get(load(ref_path))$enrichment
} else if (grepl('^snRNAseq', ref_name)) {
    results_enrichment = readRDS(ref_path)
} else {
    results_enrichment <- readRDS(ref_path)$enrichment |>
        filter(!duplicated(ensembl))
}

#   Correlate Banksy clusters with reference data

this_cor = lapply(
    t_stats,
    function(x) {
        layer_stat_cor(
            results_enrichment, modeling_results = x, model_type = "enrichment",
            top_n = 100
        )
    }
)

#   Remove 'X' from Banksy cluster names
for (i in seq_len(length(this_cor))) {
    colnames(this_cor[[i]]) = sub('^X', '', colnames(this_cor[[i]]))
}

#   Annotate clusters
annotated_clusters = lapply(
    this_cor, annotate_registered_clusters, cutoff_merge_ratio = 0.1
)

#   Make heatmaps
pdf(file.path(plot_dir, sprintf("%s.pdf", ref_name)))
for (i in seq_len(length(this_cor))) {
    if ((ref_name == 'Visium_BayesSpace_k09') && (i == 10)) {
        #   For a progress report figure, use a special ordering just for
        #   Visium BayesSpace k = 9 and Banksy res = 1
        print(
            layer_stat_cor_plot(
                this_cor[[i]], annotation = annotated_clusters[[i]],
                heatmap_legend_param = list(title = "Cor", at = c(-1, 0, 1)),
                row_order = visium_manual_row_order,
                column_order = visium_manual_col_order
            )
        )
    } else if ((ref_name == 'multiome') && (i == 10)) {
        #   For a progress report figure, use a special ordering just for
        #   multiome and Banksy res = 1
        print(
            layer_stat_cor_plot(
                this_cor[[i]], annotation = annotated_clusters[[i]],
                heatmap_legend_param = list(title = "Cor", at = c(-1, 0, 1)),
                row_order = multiome_manual_row_order
            )
        )
    } else {
        print(
            layer_stat_cor_plot(
                this_cor[[i]], annotation = annotated_clusters[[i]],
                heatmap_legend_param = list(title = "Cor", at = c(-1, 0, 1))
            )
        )
    }
}
dev.off()

saveRDS(this_cor, file = out_path)

session_info()

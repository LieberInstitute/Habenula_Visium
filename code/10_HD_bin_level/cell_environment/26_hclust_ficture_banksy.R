library(here)
library(spatialLIBD)
library(sessioninfo)
library(tidyverse)

k = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

ficture_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'registration', 'pseudobulk_spe', sprintf('%s.rds', k)
)
banksy_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'pseudobulk_spe', '1_8_cell_types.rds'
)
f_markers_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'registration', 'modeling_results', sprintf('%s.rds', k)
)
b_markers_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'modeling_results', '1_8_cell_types.rds'
)
genes_per_cluster = 70

################################################################################
#   Functions
################################################################################

read_markers = function(model_path) {
    readRDS(model_path)$enrichment |>
        as_tibble() |>
        select(ensembl, matches('^(fdr|t_stat)_')) |>
        pivot_longer(
            cols = matches('^(fdr|t_stat)_'),
            names_to = c(".value", "cluster"),
            names_pattern = "^(fdr|t_stat)_(.*)$"
        ) |>
        dplyr::rename(gene_id = ensembl) |>
        filter(fdr < 0.05, t_stat > 0) |>
        select(cluster, gene_id, fdr)
}

################################################################################
#   Main
################################################################################

#-------------------------------------------------------------------------------
#   Merge FICTURE and Banksy objects
#-------------------------------------------------------------------------------

ficture_spe = readRDS(ficture_path)
banksy_spe = readRDS(banksy_path)

common_genes = intersect(rownames(ficture_spe), rownames(banksy_spe))
ficture_spe = ficture_spe[common_genes, ]
banksy_spe = banksy_spe[common_genes, ]

colData(ficture_spe) = colData(ficture_spe)[
    , 'registration_variable', drop = FALSE
]
colData(banksy_spe) = colData(banksy_spe)[
    , 'registration_variable', drop = FALSE
]
sce = as(cbind(ficture_spe, banksy_spe), "SingleCellExperiment")
sce$sample_id = "anything"
sce$cluster = sce$registration_variable
sce$registration_variable = NULL

sce_pb = registration_pseudobulk(
    sce, var_registration = "cluster",
    var_sample_id = "sample_id", min_ncells = 1
)

#-------------------------------------------------------------------------------
#   Grab union of cluster markers
#-------------------------------------------------------------------------------

markers = rbind(read_markers(f_markers_path), read_markers(b_markers_path)) |>
    filter(gene_id %in% rownames(sce_pb)) |>
    group_by(cluster) |>
    arrange(fdr) |>
    slice_head(n = genes_per_cluster) |>
    ungroup() |>
    pull(gene_id) |>
    unique()

message(sprintf("Using %d total markers", length(markers)))

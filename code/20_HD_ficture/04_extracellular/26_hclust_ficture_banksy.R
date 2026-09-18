library(here)
library(spatialLIBD)
library(sessioninfo)
library(tidyverse)
library(stats)
library(dendextend)

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
plot_dir = here(
    'plots', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'hclust_ficture_banksy'
)
target_num_genes = 1000

dir.create(plot_dir, showWarnings = FALSE)

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

colnames(sce_pb) = colnames(sce_pb) |>
    str_replace('^anything_', '') |>
    str_replace('^X', 'Factor_')

#-------------------------------------------------------------------------------
#   Grab union of cluster markers
#-------------------------------------------------------------------------------

markers = rbind(read_markers(f_markers_path), read_markers(b_markers_path)) |>
    filter(gene_id %in% rownames(sce_pb)) |>
    group_by(cluster) |>
    arrange(fdr) |>
    #   1.3 determined empirically to approximately get target_num_genes total
    #   markers
    slice_head(n = as.integer(target_num_genes / ncol(sce_pb) * 1.3)) |>
    ungroup() |>
    pull(gene_id) |>
    unique()

message(sprintf("Using %d total markers", length(markers)))

#-------------------------------------------------------------------------------
#   Hierarchically cluster clusters transcriptionally
#-------------------------------------------------------------------------------

mat = t(logcounts(sce_pb)[markers, ])
dist_mat = as.dist(1 - cor(t(mat)))
hc = hclust(dist_mat, method = "ward.D2")

# Color labels: FICTURE clusters in blue, others in red
dend = as.dendrogram(hc)
label_colors = ifelse(startsWith(labels(dend), "Factor_"), "blue", "red")
dend = set(dend, "labels_col", label_colors)

pdf(file.path(plot_dir, sprintf("hclust_k%d.pdf", k)), width = ncol(sce_pb) / 3)
par(mar = c(12, 4, 2, 2))
plot(dend)
par(mar = c(5, 4, 4, 2))  # reset to default
dev.off()

session_info()

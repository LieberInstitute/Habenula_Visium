library(here)
library(spatialLIBD)
library(SingleCellExperiment)
library(sessioninfo)
library(tidyverse)
library(qs2)

spe_nuc_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'contamination',
    'spe', 'raw.qs2'
)
spe_cell_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'banksy',
    'leiden_res1_8.csv'
)
ct_anno_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'cluster_annotation.csv'
)
pseudo_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'contamination',
    'registration', 'sce_pb.rds'
)
model_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'contamination',
    'registration', 'modeling_results.rds'
)

dir.create(dirname(pseudo_path), showWarnings = FALSE)

################################################################################
#   Merge nuclear and cellular objects
################################################################################

spe_nuc = qs_read(spe_nuc_path)
spe_cell = readRDS(spe_cell_path)

assays(spe_cell) = list(counts = assays(spe_cell)$counts)
names(assays(spe_nuc)) = "counts"

spe_nuc = spe_nuc[
    rownames(spe_nuc) %in% rownames(spe_cell),
    colnames(spe_nuc) %in% colnames(spe_cell)
]
spe_nuc = spe_nuc[rownames(spe_cell), colnames(spe_cell)]

spe_nuc = as(spe_nuc, "SingleCellExperiment")
spe_cell = as(spe_cell, "SingleCellExperiment")

#   Another sanity check that cell IDs match between nuclear and cellular
message(
    sprintf(
        'Correlation between nuclear and cellular bin count: %.2f',
        cor(spe_nuc$bin_count, spe_cell$bin_count)
    )
)

rowData(spe_nuc) = rowData(spe_cell)
colData(spe_nuc) = colData(spe_cell)

spe_nuc$cell_id = spe_nuc$key
spe_cell$cell_id = spe_cell$key
spe_nuc$key = paste(spe_nuc$key, "nuclear", sep = "_")
spe_cell$key = paste(spe_cell$key, "cellular", sep = "_")
spe_nuc$bin_type = "nuclear"
spe_cell$bin_type = "cellular"

sce = cbind(spe_nuc, spe_cell)

#   Add in Banksy cluster
sce$banksy_cluster = tibble(key = sce$cell_id, bin_type = sce$bin_type) |>
    left_join(read_csv(cluster_path, show_col_types = FALSE), by = 'key') |>
    mutate(banksy_cluster = paste(banksy, bin_type, sep = "_")) |>
    pull(banksy_cluster)
stopifnot(!any(is.na(sce$banksy_cluster)))

################################################################################
#   Spatial registration by Banksy cluster + bin type
################################################################################

model_results = registration_wrapper(
    sce,
    var_registration = 'banksy_cluster',
    var_sample_id = 'sample_id',
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name',
    pseudobulk_rds_file = pseudo_path
)

saveRDS(model_results, model_path)

session_info()

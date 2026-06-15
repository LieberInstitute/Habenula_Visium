library(tidyverse)
library(here)
library(sessioninfo)
library(spatialLIBD)
library(rtracklayer)

spe_cell_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
spe_extra_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA',
    'extracellular', 'spe_raw.rds'
)
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'banksy',
    'leiden_res1_8.csv'
)
ct_anno_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'cluster_annotation.csv'
)
out_dir = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'cell_vs_extra_DE',
    'sce'
)
gtf_path = '/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz'
min_umi_extra = 10

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
#   Functions
################################################################################

#   Pseudobulk and clean up colData for DE later
complete_pseudobulk = function(sce) {
    sce_pb = registration_pseudobulk(
        sce, var_registration = "compartment", var_sample_id = "sample_id"
    )

    #   Preserve relevant colData and recompute certain metrics we may need
    colData(sce_pb) = colData(sce_pb)[, c('key', 'sample_id', 'compartment')]
    sce_pb$pb_sample_id = colnames(sce_pb)
    sce_pb$donor = sub('_[12]$', '', sce_pb$sample_id)
    sce_pb$donor = factor(sce_pb$donor, levels = sort(unique(sce_pb$donor)))
    sce_pb$sum_umi = unname(colSums(counts(sce_pb)))
    sce_pb$expr_chrM = colSums(
        counts(sce_pb)[which(seqnames(sce_pb) == "chrM"), , drop = FALSE]
    )
    sce_pb$expr_chrM_ratio = sce_pb$expr_chrM / sce_pb$sum_umi
    sce_pb$compartment = factor(sce_pb$compartment, levels = c('cell', 'extra'))

    return(sce_pb)
}

################################################################################
#   Load and slightly trim objects
################################################################################

spe_cell = readRDS(spe_cell_path)
spe_extra = readRDS(spe_extra_path)

assays(spe_cell) = list(counts = assays(spe_cell)$counts)
assays(spe_extra) = list(counts = assays(spe_extra)$counts)

spe_extra$key = colnames(spe_extra)
spe_extra = spe_extra[
    rowSums(assays(spe_extra)$counts > 0) > 0,
    colSums(assays(spe_extra)$counts) > min_umi_extra
]

################################################################################
#   Subset cellular and extracellular to the same genes and cells
################################################################################

shared_genes = intersect(rownames(spe_cell), rownames(spe_extra))
shared_cells = intersect(spe_cell$key, spe_extra$key)
jaccard_cells = (
    length(shared_cells) /
    length(union(spe_cell$key, spe_extra$key))
)
jaccard_genes = (
    length(shared_genes) /
    length(union(rownames(spe_cell), rownames(spe_extra)))
)
message(
    sprintf(
        'Jaccard index of cellular and extracellular cells: %.2f', jaccard_cells
    )
)
message(
    sprintf(
        'Jaccard index of cellular and extracellular genes: %.2f', jaccard_genes
    )
)

spe_cell = spe_cell[shared_genes, shared_cells]
spe_extra = spe_extra[shared_genes, shared_cells]

################################################################################
#   Merge into one minimal object, fix sequence info, and add cell types
################################################################################

spe_cell = as(spe_cell, 'SingleCellExperiment')
spe_extra = as(spe_extra, 'SingleCellExperiment')

colData(spe_cell) = colData(spe_cell)[, c('key', 'sample_id')]
colData(spe_extra) = colData(spe_cell)
spe_cell$compartment = 'cell'
spe_extra$compartment = 'extra'

sce = cbind(spe_cell, spe_extra)

#   Fix rowRanges and seqinfo, which for some reason is corrupt
gtf_genes = import(gtf_path, feature.type = "gene")
existing_rd = rowData(sce)
rowRanges(sce) = gtf_genes[match(rownames(sce), gtf_genes$gene_id)]
mcols(rowRanges(sce)) = existing_rd
rownames(sce) = rowData(sce)$gene_id

sce$cell_type = tibble(key = sce$key) |>
    left_join(read_csv(cluster_path, show_col_types = FALSE), by = 'key') |>
    dplyr::rename(cluster = banksy) |>
    left_join(read_csv(ct_anno_path, show_col_types = FALSE), by = 'cluster') |>
    pull(fine_cell_type)
stopifnot(!any(is.na(sce$cell_type)))

################################################################################
#   Pseudobulk for each cell type, including a global object, and export
################################################################################

for (this_cell_type in unique(sce$cell_type)) {
    sce_pb = complete_pseudobulk(sce[, sce$cell_type == this_cell_type])
    clean_cell_type = this_cell_type |>
        str_replace_all('\\.', '-') |>
        str_replace_all('/', '--')
    saveRDS(
        sce_pb, file.path(out_dir, sprintf('sce_%s_pb.rds', clean_cell_type))
    )
}

sce |>
    complete_pseudobulk() |>
    saveRDS(file.path(out_dir, 'sce_all_pb.rds'))

session_info()

#   The raw nuclear SPE has several problems: it's missing key spatial
#   information, has bad sample IDs, and bad rowData. Clean these up so we have
#   a generally functional nuclear SPE

library(tidyverse)
library(here)
library(rtracklayer)
library(sessioninfo)
library(SpatialExperiment)
library(qs2)

spe_in_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'contamination',
    'spe', 'raw.qs2'
)
spe_out_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'contamination',
    'spe', 'cleaned.qs2'
)
spe_cell_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
tissue_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'tissue_key_map.csv.gz'
)
gtf_path = '/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz'

################################################################################
#   Make sample ID tissue, not donor (needed for plotting and helpful for pb)
################################################################################

spe = qs_read(spe_in_path)
spe = spe[, spe$sample_id != 'H1-6FX4YN3_D1_9902']

spe$tissue_id = tibble(key = spe$key) |>
    left_join(read_csv(tissue_path, show_col_types = FALSE), by = 'key') |>
    pull(tissue_section)
spe = spe[, !is.na(spe$tissue_id)]

spe_list = list()
for (this_tissue_id in unique(spe$tissue_id)) {
    spe_list[[this_tissue_id]] = spe[, spe$tissue_id == this_tissue_id]
    spe_list[[this_tissue_id]]$sample_id = this_tissue_id
}
spe = do.call(cbind, spe_list)
spe = do.call(cbind, spe_list) # very strange bug that fixes when this is run twice
stopifnot(all(grepl('^Br', spe$sample_id)))

#-------------------------------------------------------------------------------
#   Fix other degenerate parts of the object
#-------------------------------------------------------------------------------

colnames(spatialCoords(spe)) = c("pxl_col_in_fullres", "pxl_row_in_fullres")
spatialCoords(spe)[, 'pxl_row_in_fullres'] = -1 * spatialCoords(spe)[
    , 'pxl_row_in_fullres'
]
names(assays(spe)) = "counts"
spe$exclude_overlapping = FALSE

spe_cell = readRDS(spe_cell_path)
imgData(spe) = imgData(spe_cell)[
    imgData(spe_cell)$sample_id %in% spe$sample_id,
]

spe = spe[
    rowSums(assays(spe)$counts) > 0, colSums(assays(spe)$counts) > 10
]

#   Bad rowData
gtf_genes = import(gtf_path, feature.type = "gene")
spe = spe[rownames(spe) %in% gtf_genes$gene_id, ]
rowRanges(spe) = gtf_genes[match(rownames(spe), gtf_genes$gene_id)]
rownames(spe) = rowData(spe)$gene_id

qs_save(spe, spe_out_path)

session_info()

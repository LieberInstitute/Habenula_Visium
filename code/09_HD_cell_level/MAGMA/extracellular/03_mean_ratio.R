library(SpatialExperiment)
library(tidyverse)
library(DeconvoBuddies)
library(here)
library(sessioninfo)

spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA',
    'extracellular', 'spe_norm_filtered.rds'
)
out_dir = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA',
    'extracellular', 'gene_sets'
)

dir.create(out_dir, showWarnings = FALSE)

export_set = function(spe, cell_type_col, file_tag) {
    marker_stats = get_mean_ratio(
            sce = spe, assay_name = "logcounts", cellType_col = cell_type_col,
            gene_ensembl = "gene_id", gene_name = "gene_name"
        ) |>
        filter(MeanRatio > 1) |>
        dplyr::rename(set_id = cellType.target, gene_id = gene) |>
        select(set_id, gene_id) |>
        arrange(set_id)

    message(sprintf("Marker counts for %s resolution:", file_tag))
    print(table(marker_stats$set_id))

    write_tsv(marker_stats, file.path(out_dir, sprintf("%s.tsv", file_tag)))
}

spe = readRDS(spe_path)

#   Also define mid and broad cell-type resolutions
spe$cell_type_mid = str_replace(spe$cell_type, '^([ML])Hb\\.[^_]+', '\\1Hb')
spe$cell_type_broad = str_replace(spe$cell_type, '^[ML]Hb\\.[^_]+', 'Hb')

export_set(spe, "cell_type", "fine")
export_set(spe, "cell_type_mid", "mid")
export_set(spe, "cell_type_broad", "broad")

session_info()

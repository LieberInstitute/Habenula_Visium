library(SpatialExperiment)
library(tidyverse)
library(DeconvoBuddies)
library(here)
library(sessioninfo)
library(spatialLIBD)

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

get_median_df = function(cell_type, spe, cell_type_col) {
    median_df = tibble(
        gene = rownames(spe),
        row_median = MatrixGenerics::rowMedians(
            assays(spe)$counts[, spe[[cell_type_col]] == cell_type]
        ),
        cell_type = cell_type
    )

    return(median_df)
}

median_df = bind_rows(
        lapply(
            unique(spe$cell_type), get_median_df, spe = spe,
            cell_type_col = "cell_type"
        )
    ) |>
    filter(row_median > 0)
a = findMarkers_1vAll(
        spe, assay_name = "logcounts", cellType_col = "cell_type"
    ) |>
    filter(logFC > 0, log.FDR < log(0.05)) |>
    dplyr::rename(cell_type = cellType.target) |>
    inner_join(median_df, by = c('gene', 'cell_type'))

spe_pb = registration_pseudobulk(
    spe, var_registration = "cell_type", var_sample_id = "sample_id"
)
marker_stats = get_mean_ratio(
            sce = spe_pb, assay_name = "logcounts", cellType_col = cell_type_col,
            gene_ensembl = "gene_id", gene_name = "gene_name"
        ) |>
        filter(MeanRatio > 1) |>
        dplyr::rename(set_id = cellType.target, gene_id = gene) |>
        select(set_id, gene_id) |>
        arrange(set_id)
session_info()

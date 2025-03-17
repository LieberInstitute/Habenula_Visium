library(here)
library(tidyverse)
library(SpatialExperiment)
library(HDF5Array)
library(DeconvoBuddies)
library(sessioninfo)
library(data.table)

# Match
spe_dir = here('processed-data', '09_HD_cell_level', "spe_norm_filtered")
spe = loadHDF5SummarizedExperiment(spe_dir)

out_path <- here('processed-data', '09_HD_cell_level', 'ficture_aggregate', 'markers', 'all_samples_12.rds')
dir.create(dirname(out_path)) # create the directory if it doesn't exist

#   Add ficture clusters for this value of k to the SPE
cluster = fread("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/09_HD_cell_level/ficture_aggregate/ficture_merged.csv")

cluster_col = 'FICTURE_cluster'
spe$ficture_clusters = cluster[[cluster_col]][
    match(colData(spe)$key, cluster$key)
]
spe <- spe[, !is.na(colData(spe)$ficture_clusters)]

findMarkers_1vAll(
        spe, assay_name = "logcounts", cellType_col = "ficture_clusters",
        mod = NULL
    ) |>
    dplyr::rename(t_stat = std.logFC) |>
    select(cellType.target, t_stat, gene) |>
    pivot_wider(names_from = cellType.target, values_from = t_stat) |>
    saveRDS(file = out_path)

session_info()
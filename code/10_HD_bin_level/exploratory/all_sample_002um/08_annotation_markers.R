library(here)
library(tidyverse)
library(SpatialExperiment)
library(HDF5Array)
library(DeconvoBuddies)
library(sessioninfo)
library(data.table)

# Match
spe <- readRDS("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/ficture_harmony/spe_raw.rds")
colData(spe)$barcode <- rownames(colData(spe))
colData(spe)$key<-paste0(colData(spe)$barcode,"_",colData(spe)$sample_id)
head(colData(spe))

out_path <- here('processed-data', '10_HD_bin_level', 'ficture_harmony', 'markers', 'all_samples_12.rds')
dir.create(dirname(out_path)) # create the directory if it doesn't exist

#   Add ficture clusters for this value of k to the SPE
cluster = fread("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/spe_raw/ficture/spe_cluster_all_sample.csv")

cluster_col = 'factor_K1'
spe$ficture_clusters = cluster[[cluster_col]][
    match(colData(spe)$key, cluster$key)
]
spe <- spe[, !is.na(colData(spe)$ficture_clusters)]
logcounts(spe) <- log1p(assay(spe, "counts"))

findMarkers_1vAll(
        spe, assay_name = "logcounts", cellType_col = "ficture_clusters",
        mod = NULL
    ) |>
    dplyr::rename(t_stat = std.logFC) |>
    select(cellType.target, t_stat, gene) |>
    pivot_wider(names_from = cellType.target, values_from = t_stat) |>
    saveRDS(file = out_path)

session_info()
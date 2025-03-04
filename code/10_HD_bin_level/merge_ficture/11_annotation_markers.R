library(here)
library(tidyverse)
library(SpatialExperiment)
library(HDF5Array)
library(DeconvoBuddies)
library(sessioninfo)
library(data.table)

spe_dir = here('processed-data', '10_HD_bin_level', 'spe_norm')
spe = loadHDF5SummarizedExperiment(spe_dir)

task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
i<- unique(spe$sample_id)[task_id]

out_path <- here('processed-data', '10_HD_bin_level', 'ficture', 'markers', sprintf('%s_12.rds', as.character(i)))
dir.create(dirname(out_path), showWarnings = FALSE) # create the directory if it doesn't exist

#   Add ficture clusters for this value of k to the SPE
cluster = fread(paste0("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/spe_norm_filtered/ficture/spe_cluster_",i,".csv"))

cluster_col = 'factor_K1'
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

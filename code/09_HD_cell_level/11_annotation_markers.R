library(here)
library(tidyverse)
library(SpatialExperiment)
library(HDF5Array)
library(DeconvoBuddies)
library(sessioninfo)

k = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

spe_dir = here('processed-data', '09_HD_cell_level', 'spe_norm')
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'banksy', sprintf('k%s.csv', k)
)
out_path = here(
    'processed-data', '09_HD_cell_level', 'banksy', 'markers',
    sprintf('k%s.rds', k)
)
cluster_col = 'banksy_lambda0.2'

dir.create(dirname(out_path), showWarnings = FALSE)

spe = loadHDF5SummarizedExperiment(spe_dir)

#   Add banksy clusters for this value of k to the SPE
cluster = read_csv(cluster_path, show_col_types = FALSE)
spe$banksy_clusters = cluster[[cluster_col]][
    match(as.numeric(colnames(spe)), cluster$key)
]

findMarkers_1vAll(
        spe, assay_name = "logcounts", cellType_col = "banksy_clusters",
        mod = NULL
    ) |>
    dplyr::rename(t_stat = std.logFC) |>
    mutate(cellType.target = sprintf('k%s_%s', k, cellType.target)) |>
    select(cellType.target, t_stat, gene) |>
    pivot_wider(names_from = cellType.target, values_from = t_stat) |>
    saveRDS(file = out_path)

session_info()

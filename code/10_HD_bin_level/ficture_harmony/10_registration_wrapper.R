#   Run registration_wrapper() on bin-level FICTURE clusters

library(here)
library(spatialLIBD)
library(HDF5Array)
library(sessioninfo)
library(tidyverse)
library(rtracklayer)

#   Get k from array task ID
k = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

spe_path = here(
    'processed-data', '10_HD_bin_level', 'ficture_harmony', 'spe_raw.rds'
)
cluster_path = here(
    'processed-data', '10_HD_bin_level', 'ficture_harmony',
    'bin_level_clusters.csv.gz'
)
pseudo_path = here(
    'processed-data', '10_HD_bin_level', 'ficture_harmony', 'registration',
    'pseudobulk_spe', sprintf('%s.rds', k)
)
model_path = here(
    'processed-data', '10_HD_bin_level', 'ficture_harmony', 'registration',
    'modeling_results', sprintf('%s.rds', k)
)
gtf_path = '/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf'

dir.create(dirname(pseudo_path), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(model_path), showWarnings = FALSE)

spe = readRDS(spe_path)

#   Join in FICTURE results for this k
colData(spe) = colData(spe) |>
    as_tibble() |>
    mutate(barcode = colnames(spe)) |>
    left_join(
        read_csv(cluster_path, show_col_types = FALSE),
        by = c('sample_id', 'barcode')
    ) |>
    rename(ficture = sprintf('FICTURE_k%d', k)) |>
    mutate(ficture = factor(ficture, levels = sort(unique(ficture)))) |>
    DataFrame()
spe = spe[, !is.na(spe$ficture)]

#   Add in basic rowData, which is missing due to read10xVisiumWrapper() not
#   being possible for 2um data
gtf = import(gtf_path) |>
    as.data.frame() |>
    as_tibble() |>
    filter(type == "gene") |>
    select(gene_id, gene_name)

stopifnot(all(rownames(spe) %in% gtf$gene_id))
rowData(spe) = gtf[match(rownames(spe), gtf$gene_id), ]

#   Pseudobulk
model_results = registration_wrapper(
    spe,
    var_registration = 'ficture',
    var_sample_id = 'sample_id',
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name',
    pseudobulk_rds_file = pseudo_path
)

saveRDS(model_results, model_path)

session_info()

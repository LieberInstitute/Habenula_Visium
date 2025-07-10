#   Run registration_wrapper() on bin-level FICTURE clusters
library(here)
library(spatialLIBD)
library(HDF5Array)
library(sessioninfo)
library(tidyverse)
library(rtracklayer)
library(data.table)
options(stringsAsFactors = TRUE)

#   Get k from array task ID
k = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

spe_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'ficture_harmony','spe',
    'y_clean_spe.rds'
)
cluster_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'ficture_harmony',
    'bin_level_clusters_batch.csv.gz'
)
pseudo_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'ficture_harmony',
    'registration', 'pseudobulk_spe', 'cleaning_y', sprintf('%d.rds', k)
)
model_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'ficture_harmony',
    'registration', 'modeling_results', 'cleaning_y', sprintf('%d.rds', k)
)
gtf_path = '/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf'
ficture_colnames = c('sample_id', 'barcode', sprintf('FICTURE_k%d', k))

dir.create(dirname(pseudo_path), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(model_path), showWarnings = FALSE, recursive = TRUE)

spe = readRDS(spe_path)

#   Join in FICTURE results for this k
colData(spe) = colData(spe) |>
    as_tibble() |>
    mutate(barcode = colnames(spe)) |>
    left_join(
        fread(cluster_path, select = ficture_colnames) |>
            as_tibble(),
        by = c('sample_id', 'barcode')
    ) |>
    dplyr::rename(ficture = sprintf('FICTURE_k%d', k)) |>
    mutate(ficture = factor(ficture, levels = sort(unique(ficture)))) |>
    #mutate(ficture = factor(as.character(paste0('c',ficture)), levels = sort(unique(paste0('c',ficture))))) |>
    DataFrame()

spe = spe[, !is.na(spe$ficture)]

#   Add in basic rowData, which is missing due to read10xVisiumWrapper() not
#   being possible for 2um data
gtf = import(gtf_path) |>
    as.data.frame() |>
    as_tibble() |>
    filter(type == "gene") |>
    select(gene_id, gene_name)

#   A tiny fraction genes were used in clustering by FICTURE by don't have
#   gene symbols in the GTF. Just drop these genes from spatial registration
genes_in_gtf = rownames(spe) %in% gtf$gene_id
warning(
    sprintf(
        "Dropping %d genes not in GTF (%d total)",
        sum(!genes_in_gtf), length(genes_in_gtf)
    )
)
spe = spe[genes_in_gtf, ]
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

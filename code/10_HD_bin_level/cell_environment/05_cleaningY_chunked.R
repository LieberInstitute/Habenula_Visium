library(jaffelab)
library(here)
library(tidyverse)
library(SpatialExperiment)
library(sessioninfo)
library(Matrix)

#   Get k from array task ID
k = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

spe_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'cell_environment',
    'spe_filtered.rds'
)
out_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'cell_environment',
    'cleaningY', 'temp_chunks', sprintf('%d.rds', k)
)
num_chunks_total = 1000
this_num_chunks = 20

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

spe = readRDS(spe_path)

mod = with(colData(spe), model.matrix(~ sample_id))

#   Split genes into chunks
gene_indices = split(
    seq_len(nrow(spe)), cut(seq_len(nrow(spe)), num_chunks_total,
    labels = FALSE)
)

message(Sys.time(), ' | Executing cleaningY in chunks...')
result_list = list()
for (i in (seq_len(this_num_chunks) + this_num_chunks * (k - 1))) {
    #   Grab a chunk of genes, tracking which counts begin as zeros
    counts_chunk = assays(spe)$normcounts[gene_indices[[i]], , drop = FALSE]
    zero_mask = counts_chunk == 0

    #   Regress out sample ID and forcefully preserve zeros
    cleaned = cleaningY(counts_chunk, mod, P = 1)
    cleaned[zero_mask] = 0

    #   Ensure column-sparse data format and append to list
    result_list[[i]] = as(cleaned, "CsparseMatrix")
}

message(Sys.time(), ' | Merging all chunks...')
final_result = do.call(rbind, result_list)

message(Sys.time(), ' | Saving merged chunks...')
saveRDS(final_result, out_path)

session_info()

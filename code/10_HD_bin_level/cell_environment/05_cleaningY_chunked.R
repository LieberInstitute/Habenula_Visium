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

# split the genes
gene_indices = split(
    seq_len(nrow(spe)), cut(seq_len(nrow(spe)), num_chunks_total,
    labels = FALSE)
)

result_list = list()
for (i in (seq_len(this_num_chunks) + this_num_chunks * (k - 1))) {
    counts_chunk = assays(spe)$normcounts[gene_indices[[i]], , drop = FALSE]
    zero_mask = counts_chunk == 0
    # regression
    cleaned = cleaningY(counts_chunk, mod, P = 1)
    cleaned[zero_mask] = 0
    result_list[[i]] = as(cleaned, "CsparseMatrix")
}
final_result = do.call(rbind, result_list)

saveRDS(final_result, out_path)

session_info()

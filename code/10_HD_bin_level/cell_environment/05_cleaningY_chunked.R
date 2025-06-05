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
    seq_len(nrow(spe)),
    cut(seq_len(nrow(spe)), num_chunks_total, labels = FALSE)
)

message(Sys.time(), ' | Executing cleaningY in chunks...')
result_list = list()
negative_df_list = list()
for (i in (seq_len(this_num_chunks) + this_num_chunks * (k - 1))) {
    #   Grab a chunk of genes, tracking which counts begin as zeros
    counts_chunk = assays(spe)$normcounts[gene_indices[[i]], , drop = FALSE]
    zero_mask = counts_chunk == 0

    #   Regress out sample ID, preserving initial zeros
    cleaned = cleaningY(counts_chunk, mod, P = 1)
    cleaned[zero_mask] = 0
    cleaned = as(cleaned, "CsparseMatrix")

    #   In rare cases, cleaningY can introduce negative counts. Warn about this,
    #   but shift expression to make the minimum zero for affected genes
    mins = rowMins(cleaned)
    mask = mins < 0
    if (any(mask)) {
        negative_df_list[[i]] = tibble(
            gene_symbol = rowData(spe)$symbol[gene_indices[[i]][mask]],
            min_val = mins[mask],
            mean_val = rowMeans(cleaned)[mask]
        )
        
        #   Shift counts to make the minimum zero. Conversion to dense is
        #   critical to avoid the extremely slow rowwise operation on the
        #   previously column-sparse format
        cleaned = as.matrix(cleaned)
        cleaned[mask, ] = cleaned[mask, ] + mins[mask]
        cleaned = as(cleaned, "CsparseMatrix")
    }

    result_list[[i]] = cleaned
}

#   Warn about any genes with negative counts
if (length(negative_df_list) > 0) {
    negative_df = do.call(rbind, negative_df_list)
    warning(
        "Some negative counts for genes '",
        paste(negative_df$symbol, collapse = "', '"), "'"
    )
    print(negative_df, n = nrow(negative_df))
}

message(Sys.time(), ' | Merging all chunks...')
final_result = do.call(rbind, result_list)

message(Sys.time(), ' | Saving merged chunks...')
saveRDS(final_result, out_path)

session_info()

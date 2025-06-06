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
    counts_chunk = assays(spe)$logcounts[gene_indices[[i]], , drop = FALSE]
    zero_mask = counts_chunk == 0

    #   Regress out sample ID, preserving initial zeros
    cleaned = cleaningY(counts_chunk, mod, P = 1)
    cleaned[zero_mask] = 0

    #   Remove log transform for FICTURE
    cleaned = 2 ** cleaned - 1
    cleaned = as(cleaned, "CsparseMatrix")

    mins = rowMins(cleaned)
    mask = mins < 0
    if (any(mask)) {
        #   For genes where negative counts were introduced, calculate the
        #   original minimum (nonzero) counts across bins
        counts_chunk = 2 ** as.matrix(counts_chunk)[mask, , drop = FALSE] - 1
        counts_chunk[counts_chunk == 0] = NA
        orig_mins = rowMins(counts_chunk, na.rm = TRUE)

        negative_df_list[[i]] = tibble(
            gene_symbol = rowData(spe)$symbol[gene_indices[[i]][mask]],
            min_val = mins[mask],
            orig_min_val = orig_mins,
            mean_val = rowMeans(cleaned)[mask]
        )
        
        #   For each gene, shift counts to make the minimum equal to the
        #   pre-cleaningY minimum value. Conversion to dense is critical to
        #   avoid the extremely slow rowwise operation on the previously
        #   column-sparse 'cleaned'
        cleaned = as.matrix(cleaned)
        cleaned[mask, ] = cleaned[mask, ] + (orig_mins - mins[mask])
        cleaned[as.matrix(zero_mask)] = 0
        cleaned = as(cleaned, "CsparseMatrix")
    }

    result_list[[i]] = cleaned
}

#   Warn about any genes with negative counts
if (length(negative_df_list) > 0) {
    negative_df = do.call(rbind, negative_df_list)
    warning(
        "Some negative counts for genes '",
        paste(negative_df$gene_symbol, collapse = "', '"), "'"
    )
    print(negative_df, n = nrow(negative_df))
}

message(Sys.time(), ' | Merging all chunks...')
final_result = do.call(rbind, result_list)

message(Sys.time(), ' | Saving merged chunks...')
saveRDS(final_result, out_path)

session_info()

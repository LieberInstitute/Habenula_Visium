# MAGAMA on FICTURE cluster markers (k=30)
# prepare input data

userlib <- "/users/cliu3/R/4.5"
.libPaths(c(userlib, setdiff(.libPaths(), userlib)))

find.package("dplyr")
packageVersion("dplyr")
library(dplyr)
library(duckplyr)

# Find MeanRatio marker genes
library(SpatialExperiment)
library(tidyverse)
library(DeconvoBuddies)
library(here)
library(sessioninfo)
library(spatialLIBD)

spe_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
    'spe','y_clean_spe.rds'
)
cluster_path = here(
  'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
  'bin_level_clusters_batch.parquet'
)
out_dir = here(
    'processed-data', '10_HD_bin_level', 'new_samples2','ficture_harmony', 'MAGMA', 'gene_sets'
)

mean_ratio_threshold = 1.05
max_num_genes = 200

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(num_cores) || num_cores < 1) num_cores <- 1
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
Sys.setenv(DUCKPLYR_FALLBACK_INFO = "FALSE")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

spe = readRDS(spe_path)

cluster_df = read_parquet_duckdb(cluster_path) |>
  select(sample_id, barcode, FICTURE_k30) |>
  collect()

spe$barcode = rownames(colData(spe))
spe$ficture_cluster = tibble(sample_id = spe$sample_id, barcode = spe$barcode) |>
  left_join(cluster_df, by = c("sample_id", "barcode")) |>
  pull(FICTURE_k30)

rowData(spe)$gene_id = rownames(rowData(spe))
rowData(spe)$gene_name = rowData(spe)$symbol

keep <- !is.na(colData(spe)$ficture_cluster) & colData(spe)$ficture_cluster != "NA"
spe <- spe[, keep]

export_set = function(spe, cell_type_col, file_tag) {
    #   Pseudobulking accomplishes two things:
    #       1. Circumvents expression filtering in mean ratio which is likely
    #          too aggressive for this highly sparse data
    #       2. Applies basic expression filtering instead
    spe_pb = registration_pseudobulk(
        spe, var_registration = cell_type_col, var_sample_id = "sample_id"
    )
    
    marker_stats = get_mean_ratio(
            sce = spe_pb, assay_name = "logcounts",
            cellType_col = cell_type_col, gene_ensembl = "gene_name",
            gene_name = "gene_id"
        ) |>
        filter(MeanRatio > mean_ratio_threshold) |>
        dplyr::rename(
            set_id = cellType.target,
            gene_id = gene

        ) |>
        group_by(set_id) |>
        arrange(desc(MeanRatio), .by_group = TRUE) |>
        slice_head(n = max_num_genes) |>
        ungroup() |>
        select(set_id, gene_id, gene_ensembl, MeanRatio) |>
        arrange(set_id, desc(MeanRatio))

    message(sprintf("Marker counts for %s resolution:", file_tag))
    print(table(marker_stats$set_id))

    write_tsv(marker_stats, file.path(out_dir, sprintf("%s.tsv", file_tag)))
}

# ---- export gene sets (ficture) ----
export_set(spe, "ficture_cluster",  "ficturek30")

session_info()

# X0  X1 X10 X11 X12 X13 X14 X15 X16 X17 X18 X19  X2 X20 X21 X22 X23 X24 X25 X26 
#  32  12  42   8  48  63  32  64  80  62  60  89   9  42  87  64  54  89  69  92 
# X27 X28 X29  X3  X4  X5  X6  X7  X8  X9 
#  70  76  48  38  62   6  10  21  66   8 

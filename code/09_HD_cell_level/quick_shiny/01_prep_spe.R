library(tidyverse)
library(here)
library(sessioninfo)
library(SpatialExperiment)
library(qs2)
library(scater)
library(BiocSingular)
library(BiocParallel)

spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'banksy',
    'leiden_res1_8.csv'
)
ct_anno_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'cluster_annotation.csv'
)
crawdad_in_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'crawdad', 'input_cells_k17.csv.gz'
)
svg_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'nnSVG_out',
    'merged_SVGs.txt'
)
plot_dir = here('plots', '09_HD_cell_level', 'no_secondary', 'quick_shiny')
out_dir = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'quick_shiny'
)
bad_coldata_cols = c(
    'cell_ID', 'object_id', 'labels_joint_source', 'sum_umi_capped',
    'sum_gene_capped'
)
cell_type_levels = c(
    'MHb.1', 'MHb.2', 'Excit_LHb', 'LHb.2.7', 'LHb.4', 'LHb.4/Inhib_LHb_4.2',
    'Excit.Thal/Inhib_LHb_4.2', 'Excit.Thal', 'Astrocyte', 'Endo',
    'Endo/microglia', 'Ependymal', 'Subependymal', 'Oligo', 'OPC'
)
coldata_cols = c(
    'key', 'sample_id', 'donor', 'tissue_piece', 'array_row', 'array_col',
    'bin_count', 'sum_umi', 'sum_gene', 'expr_chrM', 'expr_chrM_ratio',
    'ManualAnnotation', 'exclude_overlapping', 'sizeFactor', 'ficture_cluster',
    'banksy_cluster', 'cell_type'
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
dir.create(plot_dir, showWarnings = FALSE)
dir.create(out_dir, showWarnings = FALSE)

################################################################################
#   Join in cell types, banksy clusters, and (extracellular) FICTURE clusters
################################################################################

spe = readRDS(spe_path)

anno_df = tibble(key = spe$key) |>
    left_join(read_csv(cluster_path, show_col_types = FALSE), by = 'key') |>
    dplyr::rename(cluster = banksy) |>
    left_join(read_csv(ct_anno_path, show_col_types = FALSE), by = 'cluster') |>
    left_join(
        read_csv(crawdad_in_path, show_col_types = FALSE) |>
            select(cell_key, ficture_cluster) |>
            dplyr::rename(key = cell_key),
        by = 'key'
    )
spe$cell_type = anno_df$fine_cell_type
spe$banksy_cluster = anno_df$cluster
spe$ficture_cluster = anno_df$ficture_cluster
stopifnot(!any(is.na(spe$cell_type)))

################################################################################
#   Clean up colData
################################################################################

#   Clean up cell types
spe = spe[, spe$cell_type != 'Drop']
stopifnot(setequal(spe$cell_type, cell_type_levels))
spe$cell_type = factor(spe$cell_type, levels = cell_type_levels)

#   Remove columns we don't need
for (this_col in bad_coldata_cols) {
    spe[[this_col]] = NULL
}

#   Use factors where appropriate
spe$tissue_piece = factor(spe$tissue_piece, levels = c(1, 2))
spe$ficture_cluster = factor(
    as.integer(spe$ficture_cluster),
    levels = sort(unique(as.integer(spe$ficture_cluster)))
)
spe$banksy_cluster = factor(
    as.integer(spe$banksy_cluster),
    levels = sort(unique(as.integer(spe$banksy_cluster)))
)

#   Add donor
spe$donor = str_extract(spe$sample_id, '^Br[0-9]{4}')
spe$donor = factor(spe$donor, levels = sort(unique(spe$donor)))

#   Reorder sensibly
stopifnot(setequal(coldata_cols, colnames(colData(spe))))
colData(spe) = colData(spe)[, coldata_cols]

################################################################################
#   Add PCA and UMAP
################################################################################

svg = readLines(svg_path)
spe = runPCA(
    spe, ncomponents = 50, subset_row = svg, BSPARAM = IrlbaParam(),
    BPPARAM = MulticoreParam(num_cores)
)
spe = runUMAP(spe, subset_row = svg, BPPARAM = MulticoreParam(num_cores))

p = plotReducedDim(spe, dimred = "UMAP", colour_by = "cell_type", point_size = 1)
png(
    file.path(plot_dir, 'UMAP_cell_type.png'), width = 7, height = 7,
    units = 'in', res = 200
)
print(p)
dev.off()

################################################################################
#   Save objects
################################################################################

#   For ExperimentHub/ spatialLIBD::fetch_data()
saveRDS(spe, file.path(out_dir, 'spe_cell_habenula_atlas.rds'))

#   For the Shiny app
assays(spe) = list(logcounts = logcounts(spe))
qs_save(spe, file.path(out_dir, 'spe_shiny.qs2'), nthreads = num_cores)

session_info()

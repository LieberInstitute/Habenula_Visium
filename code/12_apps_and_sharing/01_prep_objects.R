library(tidyverse)
library(here)
library(sessioninfo)
library(SpatialExperiment)
library(spatialLIBD)
library(qs2)
library(scater)
library(BiocSingular)
library(BiocParallel)
library(rtracklayer)

spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
spe_pb_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'pseudobulk_spe', '1_8_cell_types.rds'
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
modeling_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'modeling_results', '1_8_cell_types.rds'
)
cell_map_path = here('raw-data', 'sample_info', 'hd_cell_type_map.csv')
gtf_path = '/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz'
plot_dir = here('plots', '12_apps_and_sharing', '01_prep_objects')
out_dir = here('processed-data', '12_apps_and_sharing', '01_prep_objects')
coldata_cols = c(
    'key', 'sample_id', 'donor', 'tissue_piece', 'array_row', 'array_col',
    'bin_count', 'sum_umi', 'sum_gene', 'expr_chrM', 'expr_chrM_ratio',
    'ManualAnnotation', 'exclude_overlapping', 'sizeFactor', 'ficture_cluster',
    'banksy_cluster', 'cell_type'
)
coldata_pb_cols = c(
    'key', 'sample_id', 'donor', 'tissue_piece', 'bin_count', 'ncells',
    'sum_umi', 'sum_gene', 'expr_chrM', 'expr_chrM_ratio', 'ManualAnnotation',
    'exclude_overlapping', 'cell_type'
)
sig_genes_n = 1000

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(0)

################################################################################
#   Functions
################################################################################

fix_rowRanges = function(sce, gtf_genes) {
    existing_rd = rowData(sce)
    rowRanges(sce) = gtf_genes[match(rownames(sce), gtf_genes$gene_id)]
    mcols(rowRanges(sce)) = existing_rd
    rownames(sce) = rowData(sce)$gene_id

    return(sce)
}

calc_reduced_dims = function(spe, svg, plot_path) {
    spe = runPCA(
        spe, ncomponents = 50, subset_row = svg, BSPARAM = IrlbaParam(),
        BPPARAM = MulticoreParam(num_cores)
    )
    spe = runUMAP(spe, subset_row = svg, BPPARAM = MulticoreParam(num_cores))

    p = plotReducedDim(
        spe, dimred = "UMAP", colour_by = "cell_type", point_size = 1
    )
    png(plot_path, width = 7, height = 7, units = 'in', res = 200)
    print(p)
    dev.off()
  
    return(spe)
}

################################################################################
#   Join in cell types, banksy clusters, and (extracellular) FICTURE clusters
################################################################################

spe = readRDS(spe_path)
spe_pb = readRDS(spe_pb_path)

anno_df = tibble(key = spe$key) |>
    left_join(read_csv(cluster_path, show_col_types = FALSE), by = 'key') |>
    dplyr::rename(cluster = banksy) |>
    left_join(read_csv(ct_anno_path, show_col_types = FALSE), by = 'cluster') |>
    left_join(
        read_csv(cell_map_path, show_col_types = FALSE),
        by = c('fine_cell_type' = 'old_cell_type')
    ) |>
    left_join(
        read_csv(crawdad_in_path, show_col_types = FALSE) |>
            select(cell_key, ficture_cluster) |>
            dplyr::rename(key = cell_key),
        by = 'key'
    )
spe$cell_type = anno_df$new_cell_type
spe$banksy_cluster = anno_df$cluster
spe$ficture_cluster = anno_df$ficture_cluster
stopifnot(!any(is.na(spe$cell_type[anno_df$fine_cell_type != 'Drop'])))

################################################################################
#    Fix rowRanges and seqinfo, which for some reason is corrupt
################################################################################

gtf_genes = import(gtf_path, feature.type = "gene")
spe = fix_rowRanges(spe, gtf_genes)
spe_pb = fix_rowRanges(spe_pb, gtf_genes)

################################################################################
#   Clean up colData
################################################################################

cluster_map = read_csv(cell_map_path, show_col_types = FALSE)
rename_map = stats::setNames(
    cluster_map$new_cell_type, cluster_map$old_cell_type
)

#-------------------------------------------------------------------------------
#   Cell-level SPE
#-------------------------------------------------------------------------------

cell_type_levels = read_csv(cell_map_path, show_col_types = FALSE)$new_cell_type

#   Clean up cell types
spe = spe[, !is.na(spe$cell_type)]
stopifnot(setequal(spe$cell_type, cell_type_levels))
spe$cell_type = factor(spe$cell_type, levels = cell_type_levels)

#   Add colors for the app
spe$cell_type_colors = tibble(new_cell_type = spe$cell_type) |>
    left_join(cluster_map, by = 'new_cell_type') |>
    pull(color)

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

#   Reorder (and subset) sensibly
stopifnot(all(coldata_cols %in% colnames(colData(spe))))
colData(spe) = colData(spe)[, coldata_cols]

#-------------------------------------------------------------------------------
#   Pseudobulked SPE
#-------------------------------------------------------------------------------

spe_pb$cell_type = case_when(
    spe_pb$cell_type == 'Endo.microglia' ~ 'Endo/microglia',
    spe_pb$cell_type == 'Excit.Thal.Inhib_LHb_4.2' ~ 'Excit.Thal/Inhib_LHb_4.2',
    spe_pb$cell_type == 'LHb.4.Inhib_LHb_4.2' ~ 'LHb.4/Inhib_LHb_4.2',
    TRUE ~ spe_pb$cell_type
)
spe_pb$cell_type = factor(
    dplyr::coalesce(unname(rename_map[spe_pb$cell_type]), spe_pb$cell_type),
    levels = cell_type_levels
)

#   Add colors for the app
spe_pb$cell_type_colors = tibble(new_cell_type = spe_pb$cell_type) |>
    left_join(cluster_map, by = 'new_cell_type') |>
    pull(color)

spe_pb$key = paste(spe_pb$sample_id, spe_pb$cell_type, sep = "_")

#   Recompute bin_count as sum across all cells
bin_df = tibble(
        cell_type = spe$cell_type, bin_count = spe$bin_count,
        sample_id = spe$sample_id
    ) |>
    mutate(key = paste(sample_id, cell_type, sep = "_")) |>
    group_by(key) |>
    summarize(bin_count = sum(bin_count))
spe_pb$bin_count = tibble(key = spe_pb$key) |>
    left_join(bin_df, by = 'key') |>
    pull(bin_count)

spe_pb$tissue_piece = factor(spe_pb$tissue_piece, levels = c(1, 2))

#   Add donor
spe_pb$donor = str_extract(spe_pb$sample_id, '^Br[0-9]{4}')
spe_pb$donor = factor(spe_pb$donor, levels = sort(unique(spe_pb$donor)))

#   Recompute spatialLIBD metrics
spe_pb$sum_umi = colSums(counts(spe_pb))
spe_pb$sum_gene = colSums(counts(spe_pb) > 0)
is_mito = which(seqnames(spe_pb) == 'chrM')
spe_pb$expr_chrM = colSums(counts(spe_pb)[is_mito, , drop = FALSE])
spe_pb$expr_chrM_ratio = spe_pb$expr_chrM / spe_pb$sum_umi
spe_pb$ManualAnnotation = "NA"

#   Reorder (and subset) sensibly
stopifnot(all(coldata_pb_cols %in% colnames(colData(spe_pb))))
colData(spe_pb) = colData(spe_pb)[, coldata_pb_cols]

################################################################################
#   Add PCA and UMAP
################################################################################

svg = readLines(svg_path)
spe = calc_reduced_dims(
    spe, svg, file.path(plot_dir, 'UMAP_cell_type_cell.png')
)
spe_pb = calc_reduced_dims(
    spe_pb, svg, file.path(plot_dir, 'UMAP_cell_type_pb.png')
)

################################################################################
#   Prep modeling results and 'sig_genes' for the Shiny app
################################################################################

modeling_results = readRDS(modeling_path)

#   Rename cell-type columns from old naming convention to new.
#   Columns encode "/" as "." (e.g. "Endo.microglia"), so build a lookup from
#   dot-encoded old names to new names before renaming.
suffix_rename = stats::setNames(
    unname(rename_map),
    gsub("/", ".", names(rename_map), fixed = TRUE)
)
rename_ct_cols = function(df) {
    pattern = "^(t_stat_|p_value_|fdr_|logFC_)(.+)$"
    new_names = colnames(df)
    for (i in seq_along(new_names)) {
        m = regmatches(new_names[i], regexec(pattern, new_names[i]))[[1]]
        if (length(m) == 3) {
            new_suffix = suffix_rename[m[3]]
            if (!is.na(new_suffix)) new_names[i] = paste0(m[2], new_suffix)
        }
    }
    stats::setNames(df, new_names)
}
modeling_results = lapply(modeling_results, rename_ct_cols)

spe_pb$spatialLIBD = spe_pb$cell_type
sig_genes = sig_genes_extract_all(
    n = min(sig_genes_n, nrow(spe_pb)),
    modeling_results = modeling_results, sce_layer = spe_pb
)
spe_pb$spatialLIBD = NULL

################################################################################
#   Save objects
################################################################################

sce_pb = as(spe_pb, "SingleCellExperiment")

#   For ExperimentHub/ spatialLIBD::fetch_data()
saveRDS(spe, file.path(out_dir, 'spe_cell_habenula_atlas.rds'))
saveRDS(sce_pb, file.path(out_dir, 'sce_pb_habenula_atlas.rds'))

#   For the Shiny app
assays(spe) = list(logcounts = logcounts(spe))
assays(sce_pb) = list(logcounts = logcounts(sce_pb))

qs_save(spe, file.path(out_dir, 'spe_shiny.qs2'), nthreads = num_cores)
qs_save(sce_pb, file.path(out_dir, 'sce_pb_shiny.qs2'), nthreads = num_cores)
qs_save(
    sig_genes, file.path(out_dir, 'sig_genes_shiny.qs2'), nthreads = num_cores
)
qs_save(
    modeling_results, file.path(out_dir, 'modeling_results.qs2'),
    nthreads = num_cores
)

session_info()

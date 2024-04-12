library("spatialLIBD")
library("scran")
library("tidyverse")
library("here")
library("lobstr")
library("sessioninfo")




## Create output directories
dir_rdata <- here::here("processed-data", "02_build_spe")
if (!dir.exists(dir_rdata)) { dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE) }

dir_txtdata <- here::here("processed-data", "02_build_spe/cvs_files")
if (!dir.exists(dir_txtdata)) { dir.create(dir_txtdata, showWarnings = FALSE, recursive = TRUE) }

dir_plots <- here::here("plots", "02_build_spe")
if (!dir.exists(dir_plots)) { dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE) }

## Define some info for the samples
sample_info <- data.frame(
    sample_id = c(
        "V12D07-075_C1"
    )
)
sample_info$subject <- "Br8112"
sample_info$sample_path <-
    file.path(
        here::here("processed-data", "01_spaceranger"),
        sample_info$sample_id,
        "outs"
    )
stopifnot(all(file.exists(sample_info$sample_path)))

## Define the donor info using information from
## the habenulaPilot paper
## TODO Update this info!
donor_info <- data.frame(
    subject = c("Br8112"),
    age = c(65.75),
    sex = c("F"),
    race = "EA/CAUC",
    pmi = c(31.5),
    diagnosis = c("Control"),
    rin = c(7)
)

## Combine sample info with the donor info
sample_info <- merge(sample_info, donor_info)


## Build basic SPE
Sys.time()
spe <- read10xVisiumWrapper(
    sample_info$sample_path,
    sample_info$sample_id,
    type = "sparse",
    data = "raw",
    images = c("lowres", "hires", "detected", "aligned"),
    load = TRUE#,
    #reference_gtf = here("raw-data", "genes.gtf") ## Not needed at JHPCE
)
Sys.time()
# 2024-02-22 14:43:49.823746 SpatialExperiment::read10xVisium: reading basic data from SpaceRanger
# 2024-02-22 14:43:56.531453 read10xVisiumAnalysis: reading analysis output from SpaceRanger
# 2024-02-22 14:43:56.687378 add10xVisiumAnalysis: adding analysis output from SpaceRanger
# 2024-02-22 14:43:56.86425 rtracklayer::import: reading the reference GTF file
# 2024-02-22 14:44:28.499147 adding gene information to the SPE object
# 2024-02-22 14:44:28.523581 adding information used by spatialLIBD
# [1] "2024-02-22 14:44:28 EST"

# spe@int_colData$reducedDims
# colnames(spe)
# rownames(spe)
colnames(colData(spe))

## Add the study design info
add_design <- function(spe) {
    new_col <- merge(colData(spe), sample_info)
    ## Fix order
    new_col <- new_col[match(spe$key, new_col$key), ]
    stopifnot(identical(new_col$key, spe$key))
    rownames(new_col) <- rownames(colData(spe))
    colData(spe) <-
        new_col[, -which(colnames(new_col) == "sample_path")]
    return(spe)
}
spe <- add_design(spe)

# head(colData(spe))

# ## Read in cell counts and segmentation results
# segmentations_list <-
#     lapply(sample_info$sample_id, function(sampleid) {
#         file <-
#             here(
#                 "processed-data",
#                 "spaceranger",
#                 sampleid,
#                 "outs",
#                 "spatial",
#                 "tissue_spot_counts.csv"
#             )
#         if (!file.exists(file)) {
#             return(NULL)
#         }
#         x <- read.csv(file)
#         x$key <- paste0(x$barcode, "_", sampleid)
#         return(x)
#     })
# ## Merge them (once the these files are done, this could be replaced by an rbind)
# segmentations <-
#     Reduce(function(...) {
#         merge(..., all = TRUE)
#     }, segmentations_list[lengths(segmentations_list) > 0])
#
# ## Add the information
# segmentation_match <- match(spe$key, segmentations$key)
# segmentation_info <-
#     segmentations[segmentation_match, -which(
#         colnames(segmentations) %in% c("barcode", "tissue", "row", "col", "imagerow", "imagecol", "key")
#     )]
# colData(spe) <- cbind(colData(spe), segmentation_info)


cat("Initial number of spots:", dim(spe)[2], "\n")
# Initial number of spots: 4992

## Remove genes with no data

expr <- which(rowSums(counts(spe)) > 0)
cat("Number of genes with counts:", length(expr))
#23439
no_expr <- which(rowSums(counts(spe)) == 0)
cat("Number of genes with no counts:", length(no_expr))
# [1] 13162
length(no_expr) / nrow(spe) * 100
# [1] 35.96077
spe <- spe[-no_expr, ]

## For visualizing this later with spatialLIBD
spe$overlaps_tissue <-
    factor(ifelse(spe$in_tissue, "in", "out"))

## Save with and without dropping spots outside of the tissue
spe_raw <- spe

saveRDS(spe_raw, file.path(dir_rdata, "spe_raw.rds"))

## Size in Gb
lobstr::obj_size(spe_raw)
# 206.32 MB

## Now drop the spots outside the tissue
spe <- spe_raw[, spe_raw$in_tissue]
dim(spe)
# [1] 23439  3615

cat("Spots in tissue:", dim(spe)[2], "\n")
# Initial number of spots: 3615

## Remove spots without counts
if (any(colSums(counts(spe)) == 0)) {
    message("removing spots without counts for spe")
    spe <- spe[, -which(colSums(counts(spe)) == 0)]
    dim(spe)
}


lobstr::obj_size(spe)
# 194.23 MB

saveRDS(spe, file.path(dir_rdata, "spe.rds"))


# ==============================================================================

# copied from https://github.com/LieberInstitute/spatialDLPFC/blob/14a1f253a92e43c01fec3cc3077a9b2cf9ce9fc0/code/analysis/01_build_spe/01_build_spe.R#L216-L221

## Inspect in vs outside of tissue

vis_grid_clus(
  spe = spe_raw,
  clustervar = "in_tissue",
  pdf = here(dir_plots, "all_in_tissue_grid.pdf"),
  sort_clust = FALSE,
  colors = c("TRUE" = "grey90", "FALSE" = "orange")
)

summary(spe_raw$sum_umi[which(!colData(spe_raw)$in_tissue)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 202     803    1023    1278    1499    6969

head(table(spe_raw$sum_umi[which(!colData(spe_raw)$in_tissue)]))
# 0  1  2  3  4  5
# 12 17 50 59 65 81

# 202 428 432 436 455 482 why only 1? csc
# 1   1   1   1   1   1  csc?

vis_grid_gene(
  spe = spe_raw[, which(!colData(spe_raw)$in_tissue)],
  geneid = "sum_umi",
  pdf = here::here("plots", "02_build_spe", "out_tissue_sum_umi.pdf"),
  assayname = "counts"
)

summary(spe_raw$sum_gene[which(!colData(spe_raw)$in_tissue)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 118.0   423.0   532.0   694.5   848.0  2915.0

vis_grid_gene(
  spe = spe_raw[, which(!colData(spe_raw)$in_tissue)],
  geneid = "sum_gene",
  pdf = here::here("plots", "02_build_spe", "out_tissue_sum_gene.pdf"),
  assayname = "counts"
)

summary(spe_raw$expr_chrM_ratio[which(!colData(spe_raw)$in_tissue)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.1246  0.2205  0.3227  0.3048  0.3711  0.5323

vis_grid_gene(
  spe = spe_raw[, which(!colData(spe_raw)$in_tissue)],
  geneid = "expr_chrM_ratio",
  pdf = here::here("plots", "02_build_spe", "out_tissue_expr_chrM_ratio.pdf"),
  assayname = "counts"
)

summary(spe$sum_umi)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 61    2365    3873    4340    5527   27840

vis_grid_gene(
  spe = spe,
  geneid = "sum_umi",
  pdf = here::here("plots", "02_build_spe", "in_tissue_sum_umi.pdf"),
  assayname = "counts"
)

summary(spe$sum_gene)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 56    1179    1734    1824    2290    6520

vis_grid_gene(
  spe = spe,
  geneid = "sum_gene",
  pdf = here::here("plots", "02_build_spe", "in_tissue_sum_gene.pdf"),
  assayname = "counts"
)

summary(spe$expr_chrM_ratio)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.08964 0.21649 0.24801 0.24684 0.27610 0.42706

vis_grid_gene(
  spe = spe,
  geneid = "expr_chrM_ratio",
  pdf = here::here("plots", "02_build_spe", "in_tissue_expr_chrM_ratio.pdf"),
  assayname = "counts"
)

vis_grid_gene(
  spe = spe_raw,
  geneid = "sum_umi",
  pdf = here::here("plots", "02_build_spe", "all_sum_umi.pdf"),
  assayname = "counts"
)

vis_grid_gene(
  spe = spe_raw,
  geneid = "sum_gene",
  pdf = here::here("plots", "02_build_spe", "all_sum_gene.pdf"),
  assayname = "counts"
)

vis_grid_gene(
  spe = spe_raw,
  geneid = "expr_chrM_ratio",
  pdf = here::here("plots", "02_build_spe", "all_expr_chrM_ratio.pdf"),
  assayname = "counts"
)


# ==============================================================================

## Read in the data and add additional QC metrics, followed by filtering data on edge spots
## Adapting code from https://github.com/LieberInstitute/Habenula_Visium/blob/0e020dd1a580b2bcea130057a2d7c01e5f928001/code/04_harmony_BayesSpace/01-filter_normalize.R#L29C1-L150C1

cat("Initial number of spots:", dim(spe)[2], "\n")
# Preliminary QC
spe <- spe[
  rowSums(assays(spe)$counts) > 0,
  (colSums(assays(spe)$counts) > 0) & spe$in_tissue
]

cat("Number of spots after preliminary QC:", dim(spe)[2], "\n")
## Metrics QC
metrics_qc <- function(spe) {
  qc_df <- data.frame(
    log2sum = log2(spe$sum_umi),
    log2detected = log2(spe$sum_gene),
    subsets_Mito_percent = spe$expr_chrM_ratio * 100,
    sample_id = spe$sample_id
  )

  qcfilter <- DataFrame(
    low_lib_size = isOutlier(
      qc_df$log2sum,
      type = "lower",
      log = TRUE,
      batch = qc_df$sample_id
    ),
    low_n_features = isOutlier(
      qc_df$log2detected,
      type = "lower",
      log = TRUE,
      batch = qc_df$sample_id
    ),
    high_subsets_Mito_percent = isOutlier(
      qc_df$subsets_Mito_percent,
      type = "higher",
      batch = qc_df$sample_id
    )
  )
  qcfilter$discard <-
    (qcfilter$low_lib_size |
       qcfilter$low_n_features) | qcfilter$high_subsets_Mito_percent


  spe$scran_low_lib_size_low_mito <-
    factor(
      qcfilter$low_lib_size &
        qc_df$subsets_Mito_percent < 0.5,
      levels = c("TRUE", "FALSE")
    )


  spe$scran_discard <-
    factor(qcfilter$discard, levels = c("TRUE", "FALSE"))
  spe$scran_low_lib_size <-
    factor(qcfilter$low_lib_size, levels = c("TRUE", "FALSE"))
  spe$scran_low_n_features <-
    factor(qcfilter$low_n_features, levels = c("TRUE", "FALSE"))
  spe$scran_high_subsets_Mito_percent <-
    factor(qcfilter$high_subsets_Mito_percent,
           levels = c("TRUE", "FALSE")
    )

  ## Find edge spots
  spots <- data.frame(
    row = spe$array_row,
    col = spe$array_col,
    sample_id = spe$sample_id
  )

  # edge_spots_row <-
  #   group_by(spots, sample_id, row) %>% summarize(min_col = min(col), max_col = max(col))
  edge_spots_row <-
    group_by(spots, sample_id, row) %>% mutate(min_col = min(col), max_col = max(col))
  # edge_spots_col <-
  #   group_by(spots, sample_id, col) %>% summarize(min_row = min(row), max_row = max(row))
  edge_spots_col <-
    group_by(spots, sample_id, col) %>% mutate(min_row = min(row), max_row = max(row))

  spots <-
    left_join(spots, edge_spots_row) %>% left_join(edge_spots_col)
  spots$edge_spots <-
    with(
      spots,
      row == min_row | row == max_row | col == min_col | col == max_col
    )

  head(spots, n=3)
  # row col     sample_id min_col max_col min_row max_row edge_spots
  # 1  50 102 V12D07-075_C1      34     126       6      72      FALSE
  # 2  14  94 V12D07-075_C1       0     120       6      70      FALSE
  # 3  61  97 V12D07-075_C1      43     127       5      71      FALSE

  spots$row_distance <-
    with(spots, pmin(abs(row - min_row), abs(row - max_row)))
  spots$col_distance <-
    with(spots, pmin(abs(col - min_col), abs(col - max_col)))
  ## spots$edge_distance <- with(spots, sqrt(row_distance^2 + col_distance^2))
  ## The above is from:
  ## sqrt((x_1 - x_2)^2 + (y_1 - y_2)^2)
  ## but it was wrong, here's a case the the smallest distance is on the column:
  ## sqrt(0^2 + col_distance^2) = col_distance
  spots$edge_distance <-
    with(spots, pmin(row_distance, col_distance))


  spe$edge_spots <-
    factor(spots$edge_spots, levels = c("TRUE", "FALSE"))
  spe$edge_distance <- spots$edge_distance


  spe$scran_low_lib_size_edge <-
    factor(
      qcfilter$low_lib_size &
        spots$edge_distance < 1,
      levels = c("TRUE", "FALSE")
    )

  return(spe)
}

spe <- metrics_qc(spe)
#head(spe$scran_low_lib_size_edge)

## Save object with metrics_qc()
saveRDS(spe, file.path(dir_rdata, "spe_with_scran_low_lib_size_edge.rds"))

colnames(colData(spe))

## Locate low library size spots on the edge
addmargins(table("Low_libsize_edge" = spe$scran_low_lib_size_edge))
# Low_libsize_edge
# TRUE FALSE   Sum
# 11  3604  3615

# ==============================================================================
## plot edge empty spots
## egde_distance sample plots
# edge_distance https://github.com/LieberInstitute/Visium_SPG_AD/blob/master/plots/07_spot_qc/egde_distance_wholegenome.pdf
# https://github.com/LieberInstitute/Visium_SPG_AD/blob/master/plots/07_spot_qc/scran_targeted_low_lib_size_vs_edge_distance.pdf
# https://github.com/LieberInstitute/Visium_SPG_AD/blob/master/plots/07_spot_qc/scran_targeted_low_lib_size.pdf

vis_grid_gene(
  spe = spe,
  geneid = "edge_distance",
  pdf = file.path(dir_plots, "in_tissue_egde_distance.pdf"),
  spatial = FALSE,
  point_size = 2,
  minCount = -1,
  cont_colors = viridisLite::viridis(21, direction = -1)
)


# spe_wholegenome$quality_groups <- "Pass"
# spe_wholegenome$quality_groups[spe_wholegenome$scran_discard == "TRUE"] <- "LQ: retained"
# spe_wholegenome$quality_groups[spe_wholegenome$glare] <- "LQ: glare"
# spe_wholegenome$quality_groups[spe_wholegenome$drop_low_library_edge_either] <- "LQ: low lib size & edge"
# table(spe_wholegenome$quality_groups)
# # LQ: glare LQ: low lib size & edge            LQ: retained                    Pass
# #        20                     152                     930                   37185
#
# quality_groups_colors <- c("Pass" = "grey90", "LQ: retained" = "orange", "LQ: glare" = "steelblue3", "LQ: low lib size & edge" = "violetred")
# p_list <- vis_grid_clus(
#   spe = spe_wholegenome,
#   clustervar = "quality_groups",
#   sort_clust = FALSE,
#   colors = quality_groups_colors,
#   spatial = FALSE,
#   point_size = 2,
#   return_plots = TRUE
# )
#
# pdf(file.path(dir_plots, "scran_low_lib_size_edge.pdf"), useDingbats = FALSE, height = 8 * 4, width = 9 * 3)
# print(cowplot::plot_grid(plotlist = p_list, ncol = 1, align = "hv"))
# dev.off()

# vis_grid_gene(
#   spe = spe[, which(!colData(spe)$scran_low_lib_size_edge)],
#   geneid = "edge_spots",
#   pdf = here::here("plots", "02_build_spe", "out_tissue_sum_umi_all.pdf"),
#   assayname = "counts"
# )


# ==============================================================================
## Drop spots with a low library size that are on the edge

spe <- spe[, spe$scran_low_lib_size_edge == "FALSE"]
cat(
  "Number of spots after removed low library size spots on the tissue edge:",
  dim(spe)[2],
  "\n"
)

# Number of spots after removed low library size spots on the tissue edge: 3604





# ==============================================================================

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

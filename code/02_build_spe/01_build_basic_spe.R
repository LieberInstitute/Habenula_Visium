library("spatialLIBD")
library("here")
library("lobstr")
library("sessioninfo")



## Create output directories
dir_rdata <- here::here("processed-data", "02_build_spe")
if (!dir.exists(dir_rdata)) {
    dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
}

## Define the donor info using information from
## https://github.com/LieberInstitute/Visium_SPG_AD/blob/master/raw-data/Visium_SPG_AD_ITG_MasterExcelSummarySheet.xlsx

## Sample ID
Sid <- paste0("S", rep(1:13), "_v") 
## Slide number and slide numver
array_id <- strsplit(c("C1 A1 B1 C1 D1 A1 B1 C1 D1 A1 B1 C1 D1"), "\\s+")[[1]]
sample_slide_id <- strsplit(c("V12D07-075 V13B23-285 V13B23-285 V13B23-285 V13B23-285 V13B23-281 V13B23-281 V13B23-281 V13B23-281 V14F07-340 V14F07-340 V14F07-340 V14F07-340"), "\\s+")[[1]]
sample_slide_id <- paste0(sample_slide_id, "_", array_id)
## Brain ID
brain_id <- c("Br8112", rep("Br8518", 4), rep("Br6522", 4), rep("Br9090",4)) 

## Define some info for the samples
sample_info <- data.frame(
  sample_id_short = c(Sid), # S10_v
  sample_id = c(sample_slide_id), #V14F07-340_D1
  brain_id = c(brain_id),
  age = c(65.75, rep(41.3, 8), rep(100, 4)), # need to update next week 9090
  sex = c(rep("F", 5), rep("M", 4), rep("X", 4)), # need to update next week 9090
  race = c(rep("EA/CAUC", 13)), # need to update next week 9090
  pmi = c(31.5, rep(10.5, 4), rep(30.5, 4), rep(20, 4)), # need to update next week 9090 
  diagnosis = c("Pilot", rep("Control", 12)), # need to update next week 9090
  rin = c(7, rep(6.3, 4), rep(7.4, 4), rep(5, 4)) # need to update next week 9090
)

sample_info$sample_path <-
    file.path(
        here::here("processed-data", "01_spaceranger"),
        sample_info$sample_id,
        "outs"
    )
stopifnot(all(file.exists(sample_info$sample_path)))

## Combine sample info with the donor info

sample_info[c(colnames(sample_info))] #"sample_id", "subject", "age", "sex", "race", "pmi", "diagnosis", "rin"

## Build basic SPE
Sys.time()
spe <- read10xVisiumWrapper(
    sample_info$sample_path,
    sample_info$sample_id,
    #sample_info$sample_id_short,
    type = "sparse",
    data = "raw",
    images = c("lowres", "hires", "detected", "aligned"),
    load = TRUE,
    #reference_gtf = NULL
    reference_gtf = "/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
)
Sys.time()
# 2024-02-22 14:43:49.823746 SpatialExperiment::read10xVisium: reading basic data from SpaceRanger
# 2024-02-22 14:43:56.531453 read10xVisiumAnalysis: reading analysis output from SpaceRanger
# 2024-02-22 14:43:56.687378 add10xVisiumAnalysis: adding analysis output from SpaceRanger
# 2024-02-22 14:43:56.86425 rtracklayer::import: reading the reference GTF file
# 2024-02-22 14:44:28.499147 adding gene information to the SPE object
# 2024-02-22 14:44:28.523581 adding information used by spatialLIBD
# [1] "2024-02-22 14:44:28 EST"

# class: SpatialExperiment
# dim: 36601 24960

# spe@int_colData$reducedDims
# colnames(spe)
# rownames(spe)
colnames(colData(spe))
head(spe$sample_id)
tail(spe$sample_id)
head(assays(spe)$counts[, 1:10]) # colnames=samples; rownames=genes

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

## Remove genes with no data
expr <- which(rowSums(counts(spe)) > 0)
cat("Number of genes with counts:", length(expr))

no_expr <- which(rowSums(counts(spe)) == 0)
cat("Number of genes with no counts:", length(no_expr))
# Number of genes with no counts: 13162
# Merged samples: Number of genes with no counts: 10315

cat("% genes with counts:", (length(expr) / nrow(spe)) * 100)
cat("% genes with no counts:", (length(no_expr) / nrow(spe)) * 100)


spe <- spe[-no_expr, ]

## For visualizing this later with spatialLIBD
spe$overlaps_tissue <-
    factor(ifelse(spe$in_tissue, "in", "out"))

## Save with and without dropping spots outside of the tissue
spe_raw <- spe

## Size in Gb
lobstr::obj_size(spe_raw)
# 206.32 MB
# Merged samples: 756.82 MB (5 samples)
# Merged samples: 1.21 GB (9 samples)

saveRDS(spe_raw, file.path(dir_rdata, "spe_raw.rds"))


## Now drop the spots outside the tissue
spe <- spe_raw[, spe_raw$in_tissue]
dim(spe)
cat("Spots in tissue:", dim(spe)[2], "\n")
# Spots in tissue: 3615
# Merged samples: Spots in tissue: 28023

## Remove spots without counts
if (any(colSums(counts(spe)) == 0)) {
    message("removing spots without counts for spe")
    spe <- spe[, -which(colSums(counts(spe)) == 0)]
    dim(spe)
}
dim(spe)

lobstr::obj_size(spe)
# Merged samples: 1.16 MB

saveRDS(spe, file.path(dir_rdata, "spe.rds"))





# ==============================================================================

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


# [1] "Reproducibility information:"
# > Sys.time()
# [1] "2024-04-15 12:33:49 EDT"
# > proc.time()
# user   system  elapsed
# 272.976   11.774 2942.272
# > options(width = 120)
# > session_info()
# 1.0.4       2020-04-23 [2] CRAN (R 4.3.2)
# Biobase                * 2.62.0      2023-10-24 [2] Bioconductor
# BiocFileCache            2.10.1      2023-10-26 [2] Bioconductor
# BiocGenerics           * 0.48.1      2023-11-01 [2] Bioconductor
# BiocIO                   1.12.0      2023-10-24 [2] Bioconductor
# BiocManager              1.30.22     2023-08-08 [2] CRAN (R 4.3.2)
# BiocNeighbors            1.20.2      2024-01-07 [2] Bioconductor 3.18 (R 4.3.2)
# BiocParallel             1.36.0      2023-10-24 [2] Bioconductor
# BiocSingular             1.18.0      2023-10-24 [2] Bioconductor
# BiocVersion              3.18.1      2023-11-15 [2] Bioconductor
# Biostrings               2.70.2      2024-01-28 [2] Bioconductor 3.18 (R 4.3.2)
# bit                      4.0.5       2022-11-15 [2] CRAN (R 4.3.2)
# bit64                    4.0.5       2020-08-30 [2] CRAN (R 4.3.2)
# bitops                   1.0-7       2021-04-24 [2] CRAN (R 4.3.2)
# blob                     1.2.4       2023-03-17 [2] CRAN (R 4.3.2)
# bslib                    0.6.1       2023-11-28 [2] CRAN (R 4.3.2)
# cachem                   1.0.8       2023-05-01 [2] CRAN (R 4.3.2)
# cli                      3.6.2       2023-12-11 [2] CRAN (R 4.3.2)
# codetools                0.2-19      2023-02-01 [3] CRAN (R 4.3.2)
# colorspace               2.1-0       2023-01-23 [2] CRAN (R 4.3.2)
# config                   0.3.2       2023-08-30 [2] CRAN (R 4.3.2)
# cowplot                  1.1.3       2024-01-22 [2] CRAN (R 4.3.2)
# crayon                   1.5.2       2022-09-29 [2] CRAN (R 4.3.2)
# curl                     5.2.0       2023-12-08 [2] CRAN (R 4.3.2)
# data.table               1.15.0      2024-01-30 [2] CRAN (R 4.3.2)
# DBI                      1.2.1       2024-01-12 [2] CRAN (R 4.3.2)
# dbplyr                   2.4.0       2023-10-26 [2] CRAN (R 4.3.2)
# DelayedArray             0.28.0      2023-10-24 [2] Bioconductor
# DelayedMatrixStats       1.24.0      2023-10-24 [2] Bioconductor
# digest                   0.6.34      2024-01-11 [2] CRAN (R 4.3.2)
# doParallel               1.0.17      2022-02-07 [2] CRAN (R 4.3.2)
# dotCall64                1.1-1       2023-11-28 [2] CRAN (R 4.3.2)
# dplyr                    1.1.4       2023-11-17 [2] CRAN (R 4.3.2)
# dqrng                    0.3.2       2023-11-29 [2] CRAN (R 4.3.2)
# DropletUtils             1.22.0      2023-10-24 [2] Bioconductor
# DT                       0.31        2023-12-09 [2] CRAN (R 4.3.2)
# edgeR                    4.0.14      2024-01-29 [2] Bioconductor 3.18 (R 4.3.2)
# ellipsis                 0.3.2       2021-04-29 [2] CRAN (R 4.3.2)
# ExperimentHub            2.10.0      2023-10-24 [2] Bioconductor
# fansi                    1.0.6       2023-12-08 [2] CRAN (R 4.3.2)
# fastmap                  1.1.1       2023-02-24 [2] CRAN (R 4.3.2)
# fields                   15.2        2023-08-17 [2] CRAN (R 4.3.2)
# filelock                 1.0.3       2023-12-11 [2] CRAN (R 4.3.2)
# foreach                  1.5.2       2022-02-02 [2] CRAN (R 4.3.2)
# generics                 0.1.3       2022-07-05 [2] CRAN (R 4.3.2)
# GenomeInfoDb           * 1.38.5      2023-12-28 [2] Bioconductor 3.18 (R 4.3.2)
# GenomeInfoDbData         1.2.11      2024-02-09 [2] Bioconductor
# GenomicAlignments        1.38.2      2024-01-16 [2] Bioconductor 3.18 (R 4.3.2)
# GenomicRanges          * 1.54.1      2023-10-29 [2] Bioconductor
# ggbeeswarm               0.7.2       2023-04-29 [2] CRAN (R 4.3.2)
# ggplot2                  3.4.4       2023-10-12 [2] CRAN (R 4.3.2)
# ggrepel                  0.9.5       2024-01-10 [2] CRAN (R 4.3.2)
# glue                     1.7.0       2024-01-09 [2] CRAN (R 4.3.2)
# golem                    0.4.1       2023-06-05 [2] CRAN (R 4.3.2)
# gridExtra                2.3         2017-09-09 [2] CRAN (R 4.3.2)
# gtable                   0.3.4       2023-08-21 [2] CRAN (R 4.3.2)
# HDF5Array                1.30.0      2023-10-24 [2] Bioconductor
# here                     1.0.1       2020-12-13 [2] CRAN (R 4.3.2)
# htmltools                0.5.7       2023-11-03 [2] CRAN (R 4.3.2)
# htmlwidgets              1.6.4       2023-12-06 [2] CRAN (R 4.3.2)
# httpuv                   1.6.14      2024-01-26 [2] CRAN (R 4.3.2)
# httr                     1.4.7       2023-08-15 [2] CRAN (R 4.3.2)
# interactiveDisplayBase   1.40.0      2023-10-24 [2] Bioconductor
# IRanges                * 2.36.0      2023-10-24 [2] Bioconductor
# irlba                    2.3.5.1     2022-10-03 [2] CRAN (R 4.3.2)
# iterators                1.0.14      2022-02-05 [2] CRAN (R 4.3.2)
# jquerylib                0.1.4       2021-04-26 [2] CRAN (R 4.3.2)
# jsonlite                 1.8.8       2023-12-04 [2] CRAN (R 4.3.2)
# KEGGREST                 1.42.0      2023-10-24 [2] Bioconductor
# later                    1.3.2       2023-12-06 [2] CRAN (R 4.3.2)
# lattice                  0.22-5      2023-10-24 [3] CRAN (R 4.3.2)
# lazyeval                 0.2.2       2019-03-15 [2] CRAN (R 4.3.2)
# lifecycle                1.0.4       2023-11-07 [2] CRAN (R 4.3.2)
# limma                    3.58.1      2023-10-31 [2] Bioconductor
# lobstr                   1.1.2       2022-06-22 [2] CRAN (R 4.3.2)
# locfit                   1.5-9.8     2023-06-11 [2] CRAN (R 4.3.2)
# magick                   2.8.2       2023-12-20 [2] CRAN (R 4.3.2)
# magrittr                 2.0.3       2022-03-30 [2] CRAN (R 4.3.2)
# maps                     3.4.2       2023-12-15 [2] CRAN (R 4.3.2)
# Matrix                   1.6-5       2024-01-11 [3] CRAN (R 4.3.2)
# MatrixGenerics         * 1.14.0      2023-10-24 [2] Bioconductor
# matrixStats            * 1.2.0       2023-12-11 [2] CRAN (R 4.3.2)
# memoise                  2.0.1       2021-11-26 [2] CRAN (R 4.3.2)
# mime                     0.12        2021-09-28 [2] CRAN (R 4.3.2)
# munsell                  0.5.0       2018-06-12 [2] CRAN (R 4.3.2)
# paletteer                1.6.0       2024-01-21 [2] CRAN (R 4.3.2)
# pillar                   1.9.0       2023-03-22 [2] CRAN (R 4.3.2)
# pkgconfig                2.0.3       2019-09-22 [2] CRAN (R 4.3.2)
# plotly                   4.10.4      2024-01-13 [2] CRAN (R 4.3.2)
# png                      0.1-8       2022-11-29 [2] CRAN (R 4.3.2)
# prettyunits              1.2.0       2023-09-24 [2] CRAN (R 4.3.2)
# promises                 1.2.1       2023-08-10 [2] CRAN (R 4.3.2)
# purrr                    1.0.2       2023-08-10 [2] CRAN (R 4.3.2)
# R.methodsS3              1.8.2       2022-06-13 [2] CRAN (R 4.3.2)
# R.oo                     1.26.0      2024-01-24 [2] CRAN (R 4.3.2)
# R.utils                  2.12.3      2023-11-18 [2] CRAN (R 4.3.2)
# R6                       2.5.1       2021-08-19 [2] CRAN (R 4.3.2)
# rappdirs                 0.3.3       2021-01-31 [2] CRAN (R 4.3.2)
# RColorBrewer             1.1-3       2022-04-03 [2] CRAN (R 4.3.2)
# Rcpp                     1.0.12      2024-01-09 [2] CRAN (R 4.3.2)
# RCurl                    1.98-1.14   2024-01-09 [2] CRAN (R 4.3.2)
# rematch2                 2.1.2       2020-05-01 [2] CRAN (R 4.3.2)
# restfulr                 0.0.15      2022-06-16 [2] CRAN (R 4.3.2)
# rhdf5                    2.46.1      2023-11-29 [2] Bioconductor 3.18 (R 4.3.2)
# rhdf5filters             1.14.1      2023-11-06 [2] Bioconductor
# Rhdf5lib                 1.24.2      2024-02-07 [2] Bioconductor 3.18 (R 4.3.2)
# rjson                    0.2.21      2022-01-09 [2] CRAN (R 4.3.2)
# rlang                    1.1.3       2024-01-10 [2] CRAN (R 4.3.2)
# rprojroot                2.0.4       2023-11-05 [2] CRAN (R 4.3.2)
# Rsamtools                2.18.0      2023-10-24 [2] Bioconductor
# RSQLite                  2.3.5       2024-01-21 [2] CRAN (R 4.3.2)
# rsvd                     1.0.5       2021-04-16 [2] CRAN (R 4.3.2)
# rtracklayer              1.62.0      2023-10-24 [2] Bioconductor
# S4Arrays                 1.2.0       2023-10-24 [2] Bioconductor
# S4Vectors              * 0.40.2      2023-11-23 [2] Bioconductor 3.18 (R 4.3.2)
# sass                     0.4.8       2023-12-06 [2] CRAN (R 4.3.2)
# ScaledMatrix             1.10.0      2023-10-24 [2] Bioconductor
# scales                   1.3.0       2023-11-28 [2] CRAN (R 4.3.2)
# scater                   1.30.1      2023-11-16 [2] Bioconductor
# scuttle                  1.12.0      2023-10-24 [2] Bioconductor
# sessioninfo            * 1.2.2       2021-12-06 [2] CRAN (R 4.3.2)
# shiny                  * 1.8.0       2023-11-17 [2] CRAN (R 4.3.2)
# shinyWidgets             0.8.1       2024-01-10 [2] CRAN (R 4.3.2)
# SingleCellExperiment   * 1.24.0      2023-10-24 [2] Bioconductor
# spam                     2.10-0      2023-10-23 [2] CRAN (R 4.3.2)
# SparseArray              1.2.3       2023-12-25 [2] Bioconductor 3.18 (R 4.3.2)
# sparseMatrixStats        1.14.0      2023-10-24 [2] Bioconductor
# SpatialExperiment      * 1.12.0      2023-10-24 [2] Bioconductor
# spatialLIBD            * 1.14.1      2023-11-30 [2] Bioconductor 3.18 (R 4.3.2)
# statmod                  1.5.0       2023-01-06 [2] CRAN (R 4.3.2)
# stringi                  1.8.3       2023-12-11 [2] CRAN (R 4.3.2)
# stringr                  1.5.1       2023-11-14 [2] CRAN (R 4.3.2)
# SummarizedExperiment   * 1.32.0      2023-10-24 [2] Bioconductor
# tibble                   3.2.1       2023-03-20 [2] CRAN (R 4.3.2)
# tidyr                    1.3.1       2024-01-24 [2] CRAN (R 4.3.2)
# tidyselect               1.2.0       2022-10-10 [2] CRAN (R 4.3.2)
# utf8                     1.2.4       2023-10-22 [2] CRAN (R 4.3.2)
# vctrs                    0.6.5       2023-12-01 [2] CRAN (R 4.3.2)
# vipor                    0.4.7       2023-12-18 [2] CRAN (R 4.3.2)
# viridis                  0.6.5       2024-01-29 [2] CRAN (R 4.3.2)
# viridisLite              0.4.2       2023-05-02 [2] CRAN (R 4.3.2)
# withr                    3.0.0       2024-01-16 [2] CRAN (R 4.3.2)
# XML                      3.99-0.16.1 2024-01-22 [2] CRAN (R 4.3.2)
# xtable                   1.8-4       2019-04-21 [2] CRAN (R 4.3.2)
# XVector                  0.42.0      2023-10-24 [2] Bioconductor
# yaml                     2.3.8       2023-12-11 [2] CRAN (R 4.3.2)
# zlibbioc                 1.48.0      2023-10-24 [2] Bioconductor
#
# [1] /users/csoto/R/4.3.x
# [2] /jhpce/shared/community/core/conda_R/4.3.x/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.3.x/R/lib64/R/library

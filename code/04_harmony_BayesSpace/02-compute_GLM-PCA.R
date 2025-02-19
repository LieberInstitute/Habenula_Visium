library("spatialLIBD")
library("SpatialExperiment")
library("here")
library("tidyverse")
library("scran")
library("scater")
library("scry")
library("BiocParallel")
library("BiocSingular")
library("bluster")
library("PCAtools")
library("sessioninfo")
# library("HDF5Array")

dir_rdata <- here("processed-data", "04_harmony_BayesSpace")
# filtered_in_path <- file.path(dir_rdata, "spe_qcED_spatialLIBD_log.rds")
## temporal testing
filtered_in_path <- file.path(dir_rdata, "spe_qcED_spatialLIBD_log_not_QCed.rds")
filtered_ordinary_path <- file.path(dir_rdata, "spe_qcED_spatialLIBD_log_GLM-PCA.rds") # new SPE with GLM-PCAs
filtered_hdf5_dir <- file.path(dir_rdata, "spe_qcED_spatialLIBD_log_GLM-PCA_hdf5")
dir_plots <- here("plots", "04_harmony_BayesSpace")

num_red_dims <- 50
num_cores <- 2 # Sys.getenv('SLURM_CPUS_ON_NODE')
set.seed(20240613)

## load a filtered spe object
spe <- readRDS(filtered_in_path)
#rowData(spe)

## Verified number of TRUE spots in tissue
in_tissue_spots <- sum(as.numeric(map(unique(spe$sample_id), ~ sum(spe$in_tissue[spe$sample_id == .x]))))
print(paste0(" Spots in tissue: ", in_tissue_spots))

################################################################################
#   Compute PCA
################################################################################

message(Sys.time(), " - Running modelGeneVar()")
## From
## http://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html#4_variance_modelling
dec <- modelGeneVar(spe,
    block = spe$sample_id,
    BPPARAM = MulticoreParam(num_cores)
)
colnames(dec$per.block)

## Plot gene variance in one plot for overview 

# set some initial values 
color_v <- c("red","blue","black","green","brown")
y_axis <- c(0)
x_axis <- c(0)

## Get max axis range
for (i in 1:length(colnames(dec$per.block))) {
  current <- dec$per.block[[i]]
  y_axis <- append(y_axis, max(current$total))
  x_axis <- append(x_axis, max(current$mean))
}
y_axis <- ceiling(max(y_axis))
x_axis <- ceiling(max(x_axis))

pdf(file.path(dir_plots, "scran_modelGeneVar.pdf"), useDingbats = FALSE)
plot(dec$per.block[[1]]$mean, dec$per.block[[1]]$total, 
     ylim =c(0, y_axis), xlim =c(0, x_axis),
     xlab = "Mean log-expression",
     ylab = "Variance")

for (i in 1:length(colnames(dec$per.block))) {
  current <- dec$per.block[[i]]
  curve(metadata(current)$trend(x), add=TRUE, col=color_v[i]) 
  }
legend("topright", legend = colnames(dec$per.block),
       col=c(color_v), lty=1:2, cex=0.8)
dev.off()

## Plot gene variance by sample

pdf(file.path(dir_plots, "scran_modelGeneVar_individual_plots.pdf"), useDingbats = FALSE)
mapply(function(block, blockname) {
    plot(
        block$mean,
        block$total,
        xlab = "Mean log-expression",
        ylab = "Variance",
        main = blockname
    )
    # points(metadata(block)$mean, metadata(block)$var, col="red")
    curve(metadata(block)$trend(x),
        col = "blue",
        add = TRUE
    )
}, dec$per.block, names(dec$per.block))
dev.off()

# Ordering by most interesting genes for inspection.
hvg <- mapply(function(block) {
  head(block[order(block$bio, decreasing=TRUE),], n=20) 
  }, dec$per.block)

capture.output(hvg, file = file.path(dir_rdata, "scran_Top20_hvgALL.csv"))

## Add symbol gene-ids 
hvg <- map(hvg, function(hvg_block) {
  hvg_block$gene_name <- rowData(spe)$gene_name[match(rownames(hvg_block), rownames(spe))]
  return(hvg_block)
} )

# map(hvg, head)


# get the top variable genes at different thresholds

message(Sys.time(), " - Running getTopHVGs()")
# By default getTopHVGs() retains all genes with positive values in the var.field column of stats
#     - prop define a numeric scalar specifying the proportion of genes to report as HVGs
#     _ further we subset to genes that have FDR less than or equal to fdr.threshold
top.hvgs.p1 <- getTopHVGs(dec, prop = 0.1)
print(paste("10% HVGs genes:", length(top.hvgs.p1)))
# top.hvgs.p2 <- getTopHVGs(dec, prop = 0.2)
# print(paste("Num HVGs for top 20 proportion:", length(top.hvgs.p2)))
# top.hvgs.p5 <- getTopHVGs(dec, prop = 0.5)
# print(paste("Num HVGs for top 50 proportion:", length(top.hvgs.p5)))
# 
# top.hvgs.fdr5 <- getTopHVGs(dec, fdr.threshold = 0.05)
# print(paste("Num HVGs at FDR = 0.05:", length(top.hvgs.fdr5)))
# top.hvgs.fdr1 <- getTopHVGs(dec, fdr.threshold = 0.01)
# print(paste("Num HVGs at FDR = 0.01:", length(top.hvgs.fdr1)))

save(top.hvgs.p1,file = file.path(dir_rdata, "top.hvgs.Rdata"))
# save(
#     top.hvgs.p1,
#     top.hvgs.p2,
#     top.hvgs.p5,
#     top.hvgs.fdr5,
#     top.hvgs.fdr1,
#     file = file.path(dir_rdata, "top.hvgs.Rdata")
# )

message(Sys.time(), " - Running runPCA()")
Sys.time()

# HVG by proportion of genes to report 
# 10p is our default named PCA for further analysis
spe <-
  runPCA(spe,
         subset_row = top.hvgs.p1,
         ncomponents = num_red_dims,
         name = "PCA"
  )
# spe <- # 20p
#   runPCA(spe,
#          subset_row = top.hvgs.p2,
#          ncomponents = num_red_dims,
#          name = "PCA_p2"
#   )
# spe <- # 50p
#   runPCA(spe,
#          subset_row = top.hvgs.p5,
#          ncomponents = num_red_dims,
#          name = "PCA_p5"
#   )

# HVG by Fold Discovery Rate
# spe <-
#     runPCA(spe,
#         subset_row = top.hvgs.fdr5,
#         ncomponents = num_red_dims,
#         name = "PCA_fdr5"
#     )
# spe <-
#     runPCA(spe,
#         subset_row = top.hvgs.fdr1,
#         ncomponents = num_red_dims,
#         name = "PCA_fdr1"
#     )

Sys.time()
reducedDimNames(spe)
plotReducedDim(spe, dimred = "PCA", colour_by = "sample_id") 
# plotReducedDim(spe, dimred = "PCA_p2", colour_by = "sample_id") 

# head(reducedDims(spe)$PCA_fdr1)

##   Plot all elbow plots in the same plot and add legends including hvg used and inflection point

lst_PCA_elbow <- list(PCA = length(top.hvgs.p1))
# lst_PCA_elbow <- list(
#   PCA = length(top.hvgs.p1), PCA_p2 = length(top.hvgs.p2), PCA_p5 = length(top.hvgs.p5),
#   PCA_fdr5 = length(top.hvgs.fdr5), PCA_fdr1 = length(top.hvgs.fdr1))

## Get max axis range
max_percentVar <- map(names(lst_PCA_elbow), ~ max(attr(reducedDim(spe, .x), "percentVar")))
y_axis <- ceiling(max(unlist(max_percentVar)) + 0.5)
x_axis <- num_red_dims

pdf(file.path(dir_plots, 'pca_elbow.pdf'), useDingbats = FALSE)

plot(
  attr(reducedDim(spe, names(lst_PCA_elbow[1])), "percentVar"), 
  #xlab = gsub("^PCA_", "PC_", .x), 
  xlab = "Dimension", 
  ylab = "Variance explained (%)",
  ylim =c(0, y_axis), xlim =c(0, x_axis),
  col = color_v[1],
  main = "Elbow plots")

## Build legend list for first element
percent.var <- attr(reducedDim(spe, names(lst_PCA_elbow[1])), "percentVar") 
points(percent.var, col=color_v[1]) 
#chosen.elbow <- findElbowPoint(percent.var)
hvg.threshold <- names(lst_PCA_elbow[1])
hvg.used <- paste0("(HVG = ", as.character(lst_PCA_elbow[1]), ")")
leg <- paste(hvg.threshold, hvg.used) #, ' elbow = ', chosen.elbow)
legend_label <- c(leg)

for (i in 2:length(names(lst_PCA_elbow))) {
  percent.var <- attr(reducedDim(spe, names(lst_PCA_elbow[i])), "percentVar") 
  points(percent.var, col=color_v[i]) 
  #chosen.elbow <- findElbowPoint(percent.var)
  hvg.threshold <- names(lst_PCA_elbow[i])
  hvg.used <- paste0("(HVG = ", as.character(lst_PCA_elbow[i]), ")")
  leg <- paste(hvg.threshold, hvg.used) #, ' elbow = ', chosen.elbow)
  legend_label <- append(legend_label, leg)
}
legend("topright", legend = legend_label,
       col=c(color_v), lty=1:2, cex=0.8)

dev.off()

################################################################################
#   Compute GLM-PCA
################################################################################

message(Sys.time(), " - Running devianceFeatureSelection()")
spe <- devianceFeatureSelection(spe, assay = "counts", fam = "binomial", sorted = FALSE, batch = as.factor(spe$sample_id))
spe <- devianceFeatureSelection(spe, assay = "counts", fam = "poisson", sorted = FALSE, batch = as.factor(spe$sample_id)) 
# colnames(rowData(spe))
# head(rowData(spe)$binomial_deviance)
# binomial_dev <- rowData(spe)$binomial_deviance
# length(binomial_dev[binomial_dev == 0])
# summary(binomial_dev)
# head(rowData(spe)$poisson_deviance)

## plot binomial and poison deviance in first 100 selected genes 

pdf(file.path(dir_plots, "binomial_deviance10000.pdf"))
par(mfrow = c(2,1))
p1 <- plot(
    sort(rowData(spe)$binomial_deviance, decreasing = TRUE)[1:10000],
    type = "l",
    xlab = "ranked genes",
    ylab = "binomial deviance",
    main = "Feature Selection with Binomial Deviance"
# ) + abline(v = 10, lty = 2, col = "red") + abline(v = 20, lty = 2, col = "blue") 
) + abline(v = 1000, lty = 2, col = "red") + abline(v = 2000, lty = 2, col = "blue") + abline(v = 5000, lty = 2, col = "green") 
p2 <- plot(
    sort(rowData(spe)$poisson_deviance, decreasing = TRUE)[1:10000],
    type = "l",
    xlab = "ranked genes",
    ylab = "poisson deviance",
    main = "Feature Selection with Poisson Deviance"
# ) + abline(v = 10, lty = 2, col = "red") + abline(v = 20, lty = 2, col = "blue")
) + abline(v = 1000, lty = 2, col = "red") + abline(v = 2000, lty = 2, col = "blue") + abline(v = 5000, lty = 2, col = "green") 

plts <- p1 / p2
plts
dev.off()

## calculate residuals from binomial model

message(Sys.time(), " - Running nullResiduals()")
spe <- nullResiduals( # default params
    spe,
    assay = "counts",
    fam = "binomial",
    type = "deviance"
    # batch = as.factor(spe$sample_id)
)
# produce residual vs. fitted plot. CSC 
# assayNames(spe)
# binom_dev_residuals <- assay(spe,"binomial_deviance_residuals")
# plot(binom_dev_residuals) 

## Get HVDG
hdgs.hb.1000 <-
      rownames(spe)[order(rowData(spe)$binomial_deviance, decreasing = TRUE)][1:1000]
# hdgs.hb.2000 <-
#     rownames(spe)[order(rowData(spe)$binomial_deviance, decreasing = TRUE)][1:2000]
# hdgs.hb.5000 <-
#    rownames(spe)[order(rowData(spe)$binomial_deviance, decreasing = TRUE)][1:5000]

save(hdgs.hb.1000, file = file.path(dir_rdata, "hdgs.hb.Rdata"))
# save(hdgs.hb.1000,
#     hdgs.hb.2000,
#     hdgs.hb.5000,
#     file = file.path(dir_rdata, "hdgs.hb.Rdata")
# )

message(Sys.time(), " - Running GLM-PCA")
spe <- runPCA(
    spe,
    exprs_values = "binomial_deviance_residuals",
    subset_row = hdgs.hb.1000,
    ncomponents = num_red_dims,
    name = "GLMPCA_approx",
    BSPARAM = BiocSingular::IrlbaParam()
)

# spe <- runPCA(
#     spe,
#     exprs_values = "binomial_deviance_residuals",
#     subset_row = hdgs.hb.2000,
#     ncomponents = num_red_dims,
#     name = "GLMPCA_approx_2000",
#     BSPARAM = BiocSingular::IrlbaParam()
# )

plotReducedDim(spe, dimred = "GLMPCA_approx", colour_by = "sample_id") 
# plotReducedDim(spe, dimred = "GLMPCA_approx_2000", colour_by = "sample_id") 

## Save the processed SPE object

# message(Sys.time(), " - Saving HDF5-backed filtered spe")
# spe = saveHDF5SummarizedExperiment(
#     spe, dir = filtered_hdf5_dir, replace = TRUE
# )
# spe = realize(spe)

message(Sys.time(), " - Saving ordinary filtered spe")
saveRDS(spe, filtered_ordinary_path)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


################################################################################

# > print("Reproducibility information:")
# [1] "Reproducibility information:"
# > Sys.time()
# [1] "2024-05-24 13:21:04 EDT"
# > proc.time()
# user   system  elapsed
# 1462.945   52.427 1796.425
# > options(width = 120)
# > session_info()
# [2] CRAN (R 4.3.2)
# benchmarkmeData          1.0.4       2020-04-23 [2] CRAN (R 4.3.2)
# Biobase                * 2.62.0      2023-10-24 [2] Bioconductor
# BiocFileCache            2.10.1      2023-10-26 [2] Bioconductor
# BiocGenerics           * 0.48.1      2023-11-01 [2] Bioconductor
# BiocIO                   1.12.0      2023-10-24 [2] Bioconductor
# BiocManager              1.30.22     2023-08-08 [2] CRAN (R 4.3.2)
# BiocNeighbors            1.20.2      2024-01-07 [2] Bioconductor 3.18 (R 4.3.2)
# BiocParallel           * 1.36.0      2023-10-24 [2] Bioconductor
# BiocSingular           * 1.18.0      2023-10-24 [2] Bioconductor
# BiocVersion              3.18.1      2023-11-15 [2] Bioconductor
# Biostrings               2.70.2      2024-01-28 [2] Bioconductor 3.18 (R 4.3.2)
# bit                      4.0.5       2022-11-15 [2] CRAN (R 4.3.2)
# bit64                    4.0.5       2020-08-30 [2] CRAN (R 4.3.2)
# bitops                   1.0-7       2021-04-24 [2] CRAN (R 4.3.2)
# blob                     1.2.4       2023-03-17 [2] CRAN (R 4.3.2)
# bluster                * 1.12.0      2023-10-24 [2] Bioconductor
# bslib                    0.6.1       2023-11-28 [2] CRAN (R 4.3.2)
# cachem                   1.0.8       2023-05-01 [2] CRAN (R 4.3.2)
# cli                      3.6.2       2023-12-11 [2] CRAN (R 4.3.2)
# cluster                  2.1.6       2023-12-01 [3] CRAN (R 4.3.2)
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
# dplyr                  * 1.1.4       2023-11-17 [2] CRAN (R 4.3.2)
# dqrng                    0.3.2       2023-11-29 [2] CRAN (R 4.3.2)
# DT                       0.31        2023-12-09 [2] CRAN (R 4.3.2)
# edgeR                    4.0.14      2024-01-29 [2] Bioconductor 3.18 (R 4.3.2)
# ellipsis                 0.3.2       2021-04-29 [2] CRAN (R 4.3.2)
# ExperimentHub            2.10.0      2023-10-24 [2] Bioconductor
# fansi                    1.0.6       2023-12-08 [2] CRAN (R 4.3.2)
# fastmap                  1.1.1       2023-02-24 [2] CRAN (R 4.3.2)
# fields                   15.2        2023-08-17 [2] CRAN (R 4.3.2)
# filelock                 1.0.3       2023-12-11 [2] CRAN (R 4.3.2)
# forcats                * 1.0.0       2023-01-29 [2] CRAN (R 4.3.2)
# foreach                  1.5.2       2022-02-02 [2] CRAN (R 4.3.2)
# generics                 0.1.3       2022-07-05 [2] CRAN (R 4.3.2)
# GenomeInfoDb           * 1.38.5      2023-12-28 [2] Bioconductor 3.18 (R 4.3.2)
# GenomeInfoDbData         1.2.11      2024-02-09 [2] Bioconductor
# GenomicAlignments        1.38.2      2024-01-16 [2] Bioconductor 3.18 (R 4.3.2)
# GenomicRanges          * 1.54.1      2023-10-29 [2] Bioconductor
# ggbeeswarm               0.7.2       2023-04-29 [2] CRAN (R 4.3.2)
# ggplot2                * 3.4.4       2023-10-12 [2] CRAN (R 4.3.2)
# ggrepel                  0.9.5       2024-01-10 [2] CRAN (R 4.3.2)
# glue                     1.7.0       2024-01-09 [2] CRAN (R 4.3.2)
# golem                    0.4.1       2023-06-05 [2] CRAN (R 4.3.2)
# gridExtra                2.3         2017-09-09 [2] CRAN (R 4.3.2)
# gtable                   0.3.4       2023-08-21 [2] CRAN (R 4.3.2)
# here                   * 1.0.1       2020-12-13 [2] CRAN (R 4.3.2)
# hms                      1.1.3       2023-03-21 [2] CRAN (R 4.3.2)
# htmltools                0.5.7       2023-11-03 [2] CRAN (R 4.3.2)
# htmlwidgets              1.6.4       2023-12-06 [2] CRAN (R 4.3.2)
# httpuv                   1.6.14      2024-01-26 [2] CRAN (R 4.3.2)
# httr                     1.4.7       2023-08-15 [2] CRAN (R 4.3.2)
# igraph                   2.0.1.9008  2024-02-09 [2] Github (igraph/rigraph@39158c6)
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
# locfit                   1.5-9.8     2023-06-11 [2] CRAN (R 4.3.2)
# lubridate              * 1.9.3       2023-09-27 [2] CRAN (R 4.3.2)
# magick                   2.8.2       2023-12-20 [2] CRAN (R 4.3.2)
# magrittr                 2.0.3       2022-03-30 [2] CRAN (R 4.3.2)
# maps                     3.4.2       2023-12-15 [2] CRAN (R 4.3.2)
# Matrix                   1.6-5       2024-01-11 [3] CRAN (R 4.3.2)
# MatrixGenerics         * 1.14.0      2023-10-24 [2] Bioconductor
# matrixStats            * 1.2.0       2023-12-11 [2] CRAN (R 4.3.2)
# memoise                  2.0.1       2021-11-26 [2] CRAN (R 4.3.2)
# metapod                  1.10.1      2023-12-24 [2] Bioconductor 3.18 (R 4.3.2)
# mime                     0.12        2021-09-28 [2] CRAN (R 4.3.2)
# munsell                  0.5.0       2018-06-12 [2] CRAN (R 4.3.2)
# paletteer                1.6.0       2024-01-21 [2] CRAN (R 4.3.2)
# pillar                   1.9.0       2023-03-22 [2] CRAN (R 4.3.2)
# pkgconfig                2.0.3       2019-09-22 [2] CRAN (R 4.3.2)
# plotly                   4.10.4      2024-01-13 [2] CRAN (R 4.3.2)
# png                      0.1-8       2022-11-29 [2] CRAN (R 4.3.2)
# promises                 1.2.1       2023-08-10 [2] CRAN (R 4.3.2)
# purrr                  * 1.0.2       2023-08-10 [2] CRAN (R 4.3.2)
# R6                       2.5.1       2021-08-19 [2] CRAN (R 4.3.2)
# rappdirs                 0.3.3       2021-01-31 [2] CRAN (R 4.3.2)
# RColorBrewer             1.1-3       2022-04-03 [2] CRAN (R 4.3.2)
# Rcpp                     1.0.12      2024-01-09 [2] CRAN (R 4.3.2)
# RCurl                    1.98-1.14   2024-01-09 [2] CRAN (R 4.3.2)
# readr                  * 2.1.5       2024-01-10 [2] CRAN (R 4.3.2)
# rematch2                 2.1.2       2020-05-01 [2] CRAN (R 4.3.2)
# restfulr                 0.0.15      2022-06-16 [2] CRAN (R 4.3.2)
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
# scater                 * 1.30.1      2023-11-16 [2] Bioconductor
# scran                  * 1.30.2      2024-01-22 [2] Bioconductor 3.18 (R 4.3.2)
# scry                   * 1.14.0      2023-10-24 [2] Bioconductor
# scuttle                * 1.12.0      2023-10-24 [2] Bioconductor
# sessioninfo            * 1.2.2       2021-12-06 [2] CRAN (R 4.3.2)
# shiny                    1.8.0       2023-11-17 [2] CRAN (R 4.3.2)
# shinyWidgets             0.8.1       2024-01-10 [2] CRAN (R 4.3.2)
# SingleCellExperiment   * 1.24.0      2023-10-24 [2] Bioconductor
# spam                     2.10-0      2023-10-23 [2] CRAN (R 4.3.2)
# SparseArray              1.2.3       2023-12-25 [2] Bioconductor 3.18 (R 4.3.2)
# sparseMatrixStats        1.14.0      2023-10-24 [2] Bioconductor
# SpatialExperiment      * 1.12.0      2023-10-24 [2] Bioconductor
# spatialLIBD            * 1.15.4      2024-05-01 [1] Github (LieberInstitute/spatialLIBD@77a5303)
# statmod                  1.5.0       2023-01-06 [2] CRAN (R 4.3.2)
# stringi                  1.8.3       2023-12-11 [2] CRAN (R 4.3.2)
# stringr                * 1.5.1       2023-11-14 [2] CRAN (R 4.3.2)
# SummarizedExperiment   * 1.32.0      2023-10-24 [2] Bioconductor
# tibble                 * 3.2.1       2023-03-20 [2] CRAN (R 4.3.2)
# tidyr                  * 1.3.1       2024-01-24 [2] CRAN (R 4.3.2)
# tidyselect               1.2.0       2022-10-10 [2] CRAN (R 4.3.2)
# tidyverse              * 2.0.0       2023-02-22 [2] CRAN (R 4.3.2)
# timechange               0.3.0       2024-01-18 [2] CRAN (R 4.3.2)
# tzdb                     0.4.0       2023-05-12 [2] CRAN (R 4.3.2)
# utf8                     1.2.4       2023-10-22 [2] CRAN (R 4.3.2)
# vctrs                    0.6.5       2023-12-01 [2] CRAN (R 4.3.2)
# vipor                    0.4.7       2023-12-18 [2] CRAN (R 4.3.2)

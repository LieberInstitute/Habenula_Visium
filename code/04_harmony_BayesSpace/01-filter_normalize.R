library("spatialLIBD")
library("here")
library("tidyverse")
library("scran")
library("BiocParallel")
library("scater")
library("scry")
library("BiocSingular")
library("sessioninfo")
# library("HDF5Array")

dir_rdata <- here("processed-data", "04_harmony_BayesSpace")
# raw_in_path <- here("processed-data", "02_build_spe", "spe.rds")
filtered_in_path <- here("processed-data", "02_build_spe", "spe_qc_low_lib_edge.rds") # this is a filtered spe object
filtered_ordinary_path <- file.path(dir_rdata, "spe_filtered.rds")
filtered_hdf5_dir <- file.path(dir_rdata, "spe_filtered_hdf5")
dir_plots <- here("plots", "04_harmony_BayesSpace")
num_red_dims <- 50

num_cores <- 2 # Sys.getenv('SLURM_CPUS_ON_NODE')
set.seed(20240223)

## Create output directories
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

## load a filtered spe object
spe <- readRDS(filtered_in_path)
## Verified number of TRUE spots in tissue
in_tissue_spots <- sum(as.numeric(map(unique(spe$sample_id), ~ sum( spe$in_tissue[spe$sample_id == .x]))))
print(paste0(' Spots in tissue: ', in_tissue_spots))


################################################################################
#   Read in the data and add additional QC metrics, followed by filtering data
################################################################################
# spe <- readRDS(raw_in_path)
# cat("Initial number of spots:", dim(spe)[2], "\n")
# # Preliminary QC
# spe <- spe[
#     rowSums(assays(spe)$counts) > 0,
#     (colSums(assays(spe)$counts) > 0) & spe$in_tissue
# ]
# 
# cat("Number of spots after preliminary QC:", dim(spe)[2], "\n")
# ## Metrics QC
# metrics_qc <- function(spe) {
#     qc_df <- data.frame(
#         log2sum = log2(spe$sum_umi),
#         log2detected = log2(spe$sum_gene),
#         subsets_Mito_percent = spe$expr_chrM_ratio * 100,
#         sample_id = spe$sample_id
#     )
# 
#     qcfilter <- DataFrame(
#         low_lib_size = isOutlier(
#             qc_df$log2sum,
#             type = "lower",
#             log = TRUE,
#             batch = qc_df$sample_id
#         ),
#         low_n_features = isOutlier(
#             qc_df$log2detected,
#             type = "lower",
#             log = TRUE,
#             batch = qc_df$sample_id
#         ),
#         high_subsets_Mito_percent = isOutlier(
#             qc_df$subsets_Mito_percent,
#             type = "higher",
#             batch = qc_df$sample_id
#         )
#     )
#     qcfilter$discard <-
#         (qcfilter$low_lib_size |
#             qcfilter$low_n_features) | qcfilter$high_subsets_Mito_percent
# 
# 
#     spe$scran_low_lib_size_low_mito <-
#         factor(
#             qcfilter$low_lib_size &
#                 qc_df$subsets_Mito_percent < 0.5,
#             levels = c("TRUE", "FALSE")
#         )
# 
# 
#     spe$scran_discard <-
#         factor(qcfilter$discard, levels = c("TRUE", "FALSE"))
#     spe$scran_low_lib_size <-
#         factor(qcfilter$low_lib_size, levels = c("TRUE", "FALSE"))
#     spe$scran_low_n_features <-
#         factor(qcfilter$low_n_features, levels = c("TRUE", "FALSE"))
#     spe$scran_high_subsets_Mito_percent <-
#         factor(qcfilter$high_subsets_Mito_percent,
#             levels = c("TRUE", "FALSE")
#         )
# 
#     ## Find edge spots
#     spots <- data.frame(
#         row = spe$array_row,
#         col = spe$array_col,
#         sample_id = spe$sample_id
#     )
# 
#     edge_spots_row <-
#         group_by(spots, sample_id, row) %>% summarize(min_col = min(col), max_col = max(col))
#     edge_spots_col <-
#         group_by(spots, sample_id, col) %>% summarize(min_row = min(row), max_row = max(row))
# 
#     spots <-
#         left_join(spots, edge_spots_row) %>% left_join(edge_spots_col)
#     spots$edge_spots <-
#         with(
#             spots,
#             row == min_row | row == max_row | col == min_col | col == max_col
#         )
# 
#     spots$row_distance <-
#         with(spots, pmin(abs(row - min_row), abs(row - max_row)))
#     spots$col_distance <-
#         with(spots, pmin(abs(col - min_col), abs(col - max_col)))
#     ## spots$edge_distance <- with(spots, sqrt(row_distance^2 + col_distance^2))
#     ## The above is from:
#     ## sqrt((x_1 - x_2)^2 + (y_1 - y_2)^2)
#     ## but it was wrong, here's a case the the smallest distance is on the column:
#     ## sqrt(0^2 + col_distance^2) = col_distance
#     spots$edge_distance <-
#         with(spots, pmin(row_distance, col_distance))
# 
# 
#     spe$edge_spots <-
#         factor(spots$edge_spots, levels = c("TRUE", "FALSE"))
#     spe$edge_distance <- spots$edge_distance
# 
# 
#     spe$scran_low_lib_size_edge <-
#         factor(
#             qcfilter$low_lib_size &
#                 spots$edge_distance < 1,
#             levels = c("TRUE", "FALSE")
#         )
# 
#     return(spe)
# }
# 
# spe <- metrics_qc(spe)
# 
# ## Save object with metrics_qc()
# # saveRDS(spe, file.path(dir_rdata, "spe_with_scran_low_lib_size_edge.rds"))
# 
# ## Drop spots with a low library size that are on the edge
# spe <- spe[, spe$scran_low_lib_size_edge == "FALSE"]
# cat(
#     "Number of spots after removed low library size spots on the tissue edge:",
#     dim(spe)[2],
#     "\n"
# )

################################################################################
#   Compute log-normalized counts
################################################################################



#   Filter SPE: take only spots in tissue, drop spots with 0 counts for all
#   genes, and drop genes with 0 counts in every spot
message(Sys.time(), " - Running quickCluster()")

Sys.time()
spe$scran_quick_cluster <- quickCluster(
    spe,
    BPPARAM = MulticoreParam(num_cores),
    block = spe$sample_id,
    block.BPPARAM = MulticoreParam(num_cores)
)
Sys.time()

message(Sys.time(), " - Running computeSumFactors()")
Sys.time()
spe <- computeSumFactors(spe,
    # clusters = spe$scran_quick_cluster,
    BPPARAM = MulticoreParam(num_cores)
)
Sys.time()

print("Quick cluster table:")
table(spe$scran_quick_cluster)

message(Sys.time(), " - Running checking sizeFactors()")
summary(sizeFactors(spe))

message(Sys.time(), " - Running logNormCounts()")
spe <- logNormCounts(spe)

# #   Save a copy of the SPE with HDF5-backed assays, which will be important to
# #   control memory consumption later
# message(Sys.time(), " - Saving HDF5-backed object to control memory later")
# spe = saveHDF5SummarizedExperiment(
#     spe, dir = paste0(filtered_hdf5_dir, '_temp'), replace = TRUE
# )
# gc()

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

pdf(file.path(dir_plots, "scran_modelGeneVar.pdf"), useDingbats = FALSE)
mapply(function(block, blockname) {
    plot(
        block$mean,
        block$total,
        xlab = "Mean log-expression",
        ylab = "Variance",
        main = blockname
    )
    curve(metadata(block)$trend(x),
        col = "blue",
        add = TRUE
    )
}, dec$per.block, names(dec$per.block))
dev.off()

message(Sys.time(), " - Running getTopHVGs()")
top.hvgs.p1 <- getTopHVGs(dec, prop = 0.1)
print(paste("Num HVGs for top 10% prop:", length(top.hvgs.p1)))
top.hvgs.p2 <- getTopHVGs(dec, prop = 0.2)
print(paste("Num HVGs for top 20% prop:", length(top.hvgs.p2)))
top.hvgs.p5 <- getTopHVGs(dec, prop = 0.5)
print(paste("Num HVGs for top 50% prop:", length(top.hvgs.p5)))

top.hvgs.fdr5 <- getTopHVGs(dec, fdr.threshold = 0.05)
print(paste("Num HVGs at FDR = 0.05:", length(top.hvgs.fdr5)))

top.hvgs.fdr1 <- getTopHVGs(dec, fdr.threshold = 0.01)
print(paste("Num HVGs at FDR = 0.01:", length(top.hvgs.fdr1)))

save(
    top.hvgs.p1,
    top.hvgs.p2,
    top.hvgs.p5,
    top.hvgs.fdr5,
    top.hvgs.fdr1,
    file = file.path(dir_rdata, "top.hvgs.Rdata")
)

message(Sys.time(), " - Running runPCA()")
Sys.time()
spe <-
    runPCA(spe,
        subset_row = top.hvgs.fdr5,
        ncomponents = num_red_dims,
        name = "PCA"
    )
spe <-
    runPCA(spe,
        subset_row = top.hvgs.p1,
        ncomponents = num_red_dims,
        name = "PCA_p1"
    )
spe <-
    runPCA(spe,
        subset_row = top.hvgs.p2,
        ncomponents = num_red_dims,
        name = "PCA_p2"
    )
spe <-
    runPCA(spe,
        subset_row = top.hvgs.p5,
        ncomponents = num_red_dims,
        name = "PCA_p5"
    )
Sys.time()

#   Plot variance explained
percent.var <- attr(reducedDim(spe, "PCA_p1"), "percentVar")
pdf(file.path(dir_plots, "pca_elbow_p1.pdf"), useDingbats = FALSE)
plot(percent.var, xlab = "PC_p1", ylab = "Variance explained (%)")
dev.off()

percent.var <- attr(reducedDim(spe, "PCA_p2"), "percentVar")
pdf(file.path(dir_plots, "pca_elbow_p2.pdf"), useDingbats = FALSE)
plot(percent.var, xlab = "PC_p2", ylab = "Variance explained (%)")
dev.off()

percent.var <- attr(reducedDim(spe, "PCA_p5"), "percentVar")
pdf(file.path(dir_plots, "pca_elbow_p5.pdf"), useDingbats = FALSE)
plot(percent.var, xlab = "PC_p5", ylab = "Variance explained (%)")
dev.off()


################################################################################
#   Compute GLM-PCA
################################################################################

message(Sys.time(), " - Running devianceFeatureSelection()")
spe <- devianceFeatureSelection(
    spe,
    assay = "counts",
    fam = "binomial",
    sorted = FALSE,
    batch = as.factor(spe$sample_id)
)

pdf(file.path(dir_plots, "binomial_deviance.pdf"))
plot(
    sort(rowData(spe)$binomial_deviance, decreasing = T),
    type = "l",
    xlab = "ranked genes",
    ylab = "binomial deviance",
    main = "Feature Selection with Deviance"
)
abline(v = 2000, lty = 2, col = "red")
dev.off()

message(Sys.time(), " - Running nullResiduals()")
spe <- nullResiduals( # default params
    spe,
    assay = "counts",
    fam = "binomial",
    type = "deviance"
)


hdgs.hb.2000 <-
    rownames(spe)[order(rowData(spe)$binomial_deviance, decreasing = TRUE)][1:2000]
hdgs.hb.5000 <-
    rownames(spe)[order(rowData(spe)$binomial_deviance, decreasing = TRUE)][1:5000]
hdgs.hb.10000 <-
    rownames(spe)[order(rowData(spe)$binomial_deviance, decreasing = TRUE)][1:10000]

save(hdgs.hb.2000,
    hdgs.hb.5000,
    hdgs.hb.10000,
    file = file.path(dir_rdata, "hdgs.hb.Rdata")
)

message(Sys.time(), " - Running GLM-PCA")
spe <- runPCA(
    spe,
    exprs_values = "binomial_deviance_residuals",
    subset_row = hdgs.hb.2000,
    ncomponents = num_red_dims,
    name = "GLMPCA_approx",
    BSPARAM = BiocSingular::IrlbaParam()
)

spe <- runPCA(
    spe,
    exprs_values = "binomial_deviance_residuals",
    subset_row = hdgs.hb.5000,
    ncomponents = num_red_dims,
    name = "GLMPCA_approx_5000",
    BSPARAM = BiocSingular::IrlbaParam()
)

spe <- runPCA(
    spe,
    exprs_values = "binomial_deviance_residuals",
    subset_row = hdgs.hb.10000,
    ncomponents = num_red_dims,
    name = "GLMPCA_approx_10000",
    BSPARAM = BiocSingular::IrlbaParam()
)

################################################################################
#   Obtain preliminary clusters based on default GLM-PCA and PCA settings
################################################################################
spe$leiden20_PCA <- clusterCells(spe,
    use.dimred = "PCA",
    BLUSPARAM = SNNGraphParam(
        k = 20,
        cluster.fun = "leiden"
    )
)

spe$leiden20_GLMPCA <- clusterCells(spe,
    use.dimred = "GLMPCA_approx",
    BLUSPARAM = SNNGraphParam(
        k = 20,
        cluster.fun = "leiden"
    )
)


################################################################################
#   Save the processed SPE object
################################################################################

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

library("spatialLIBD")
library("here")
library("tidyverse")
library("scran")
library("BiocParallel")
library("scater")
library("scry")
library("BiocSingular")
library("sessioninfo")
library("HDF5Array")

dir_rdata <- here("processed-data", "04_harmony_BayesSpace")
filtered_in_path <- file.path(dir_rdata, "spe_qc_filtered_logcounts.rds")
#filtered_hdf5_dir <- file.path(dir_rdata, "spe_filtered_hdf5")
dir_plots <- here("plots", "04_harmony_BayesSpace")

num_red_dims <- 50
num_cores <- 2 # Sys.getenv('SLURM_CPUS_ON_NODE')
set.seed(20240223)

# ## Create output directories
# dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
# dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

## load a filtered spe object
spe <- readRDS(filtered_in_path)

## Verified number of TRUE spots in tissue
in_tissue_spots <- sum(as.numeric(map(unique(spe$sample_id), ~ sum( spe$in_tissue[spe$sample_id == .x]))))
print(paste0(' Spots in tissue: ', in_tissue_spots))


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

#plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
#curve(metadata(dec)$trend(x), col="blue", add=TRUE)

pdf(file.path(dir_plots, "scran_modelGeneVar.pdf"), useDingbats = FALSE)
mapply(function(block, blockname) {
    plot(
        block$mean,
        block$total,
        xlab = "Mean log-expression",
        ylab = "Variance",
        main = blockname
    )
    #points(metadata(block)$mean, metadata(block)$var, col="red")
    curve(metadata(block)$trend(x),
        col = "blue",
        add = TRUE
    )
}, dec$per.block, names(dec$per.block))
dev.off()

message(Sys.time(), " - Running getTopHVGs()")
# By default getTopHVGs() retains all genes with positive values in the var.field column of stats
#     - prop define a numeric scalar specifying the proportion of genes to report as HVGs
#     _ further we subset to genes that have FDR less than or equal to fdr.threshold
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

# lst_top_hvgs <- list(PCA_fdr1 = as.character(top.hvgs.fdr1[!is.na(top.hvgs.fdr1)]),
#                      PCA_fdr5 = as.character(top.hvgs.fdr5[!is.na(top.hvgs.fdr5)]))
# typeof(top.hvgs.fdr1)
# anyNA(top.hvgs.fdr1)
# top.hvgs.fdr1[1][1:4]
# top.hvgs.fdr1[2][1:4]
# 
# for (x in 1:length(lst_top_hvgs)) {
#   message(' Processing ', names(lst_top_hvgs[x]), lst_top_hvgs[x][1:10])
#   # spe <- runPCA(spe,
#   #     subset_row = lst_top_hvgs[x],
#   #     ncomponents = num_red_dims,
#   #     name = names(lst_top_hvgs[x])
#   # )
# }


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

# lst_PCA_elbow <- list(PCA_p1 = 'pca_elbow_p1.pdf', PCA_p2 = 'pca_elbow_p2.pdf', PCA_p5 = 'pca_elbow_p5.pdf')
# plt_elbow <- function(dim_n, plt_name) {
#   percent.var <- attr(reducedDim(spe, dim_n), "percentVar")
#   pdf(file.path(dir_plots, plt_name), useDingbats = FALSE)
#   plot(percent.var, xlab = gsub("^PCA_", "PC_", dim_n), ylab = "Variance explained (%)")
#   dev.off()
# }
# 
# map2(names(lst_PCA_elbow), lst_PCA_elbow, ~ plt_elbow(.x, .y))

percent.var <- attr(reducedDim(spe, "PCA_p1"), "percentVar")
pdf(file.path(dir_plots, "pca_elbow_p1.pdf"), useDingbats = FALSE)
plot(percent.var, xlab = "PC_p1", ylab = "Variance explained (%)")
dev.off()
#
percent.var <- attr(reducedDim(spe, "PCA_p2"), "percentVar")
pdf(file.path(dir_plots, "pca_elbow_p2.pdf"), useDingbats = FALSE)
plot(percent.var, xlab = "PC_p2", ylab = "Variance explained (%)")
dev.off()
#
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
    file = file.path(dir_rdata, "hdgs.hb_2.Rdata")
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

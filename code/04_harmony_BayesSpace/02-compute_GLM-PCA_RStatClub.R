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

here::here()

dir_rdata <- here("processed-data", "04_harmony_BayesSpace")
filtered_in_path <- file.path(dir_rdata, "spe_qcED_spatialLIBD_log.rds")
# filtered_ordinary_path <- file.path(dir_rdata, "spe_qcED_spatialLIBD_log_GLM-PCA.rds") # new SPE with GLM-PCAs
# filtered_hdf5_dir <- file.path(dir_rdata, "spe_qcED_spatialLIBD_log_GLM-PCA_hdf5")
# dir_plots <- here("plots", "04_harmony_BayesSpace")

num_red_dims <- 50
num_cores <- 2 # Sys.getenv('SLURM_CPUS_ON_NODE')
set.seed(20240613)

## load a filtered spe object
spe <- readRDS(filtered_in_path)
spe

## Verified number of TRUE spots in tissue
in_tissue_spots <- sum(as.numeric(map(unique(spe$sample_id), ~ sum(spe$in_tissue[spe$sample_id == .x]))))
print(paste0(" Spots in tissue: ", in_tissue_spots))

################################################################################
#   Compute PCA
################################################################################

if (b_hvg) {

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
  # dev.off()
  
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
  # dev.off()
  
  # Ordering by most interesting genes for inspection.
  hvg <- mapply(function(block) {
    head(block[order(block$bio, decreasing=TRUE),], n=20) 
    }, dec$per.block)
  
  capture.output(hvg, file = file.path(dir_rdata, "scran_Top20_hvgALL.csv"))
  
  # get the top variable genes at different thresholds
  
  message(Sys.time(), " - Running getTopHVGs()")
  # By default getTopHVGs() retains all genes with positive values in the var.field column of stats
  #     - prop define a numeric scalar specifying the proportion of genes to report as HVGs
  #     _ further we subset to genes that have FDR less than or equal to fdr.threshold
  top.hvgs.p1 <- getTopHVGs(dec, prop = 0.1)
  print(paste("Num HVGs for top 10 proportion:", length(top.hvgs.p1)))
  top.hvgs.p2 <- getTopHVGs(dec, prop = 0.2)
  
  # save(top.hvgs.p1, top.hvgs.p2, file = file.path(dir_rdata, "top.hvgs.Rdata"))
  
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
  spe <- # 20p
    runPCA(spe,
           subset_row = top.hvgs.p2,
           ncomponents = num_red_dims,
           name = "PCA_p2"
    )
  Sys.time()
  reducedDimNames(spe)
  plotReducedDim(spe, dimred = "PCA", colour_by = "sample_id") 
  plotReducedDim(spe, dimred = "PCA_p2", colour_by = "sample_id") 
  
  # head(reducedDims(spe)$PCA_fdr1)
  
  ##   Plot all elbow plots in the same plot and add legends including hvg used and inflection point
  
  lst_PCA_elbow <- list(
    PCA = length(top.hvgs.p1), PCA_p2 = length(top.hvgs.p2), PCA_p5 = length(top.hvgs.p5),
    PCA_fdr5 = length(top.hvgs.fdr5), PCA_fdr1 = length(top.hvgs.fdr1))
  
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
  # dev.off()

}



################################################################################
#   Compute GLM-PCA from scry
################################################################################

message(Sys.time(), " - Running devianceFeatureSelection()")

## devianceFeatureSelection(): Computes a deviance statistic for each gene for count data based on a multinomial null model that assumes each feature has a constant rate. Features with large deviance are likely to be informative. Uninformative, low deviance features can be discarded to speed up downstream analyses and reduce memory footprint.


# First we will rank genes based on deviance, to help identify the most biologically informative genes. The actual deviance values are stored in the rowData of the SingleCellExperiment object

spe <- devianceFeatureSelection(spe,                    
                                assay = "counts", 
                                fam = "binomial", 
                                sorted = FALSE, 
                                batch = as.factor(spe$sample_id))

## It returns 
# The new column name will be either binomial_deviance or poisson_deviance. 
# If the input was a matrix-like object, output is a numeric vector containing the deviance statistics for each row.

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
# dev.off()

# We can see that the deviance drops sharply after about 2,000 genes. The remaining genes are probably not informative so we discard them to speed up downstream analysis.

sce2<-sce[1:1000, ]

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
hdgs.hb.2000 <-
    rownames(spe)[order(rowData(spe)$binomial_deviance, decreasing = TRUE)][1:2000]

# save(hdgs.hb.1000, hdgs.hb.2000, file = file.path(dir_rdata, "hdgs.hb.Rdata"))

message(Sys.time(), " - Running GLM-PCA")
spe <- runPCA(
    spe,
    exprs_values = "binomial_deviance_residuals",
    subset_row = hdgs.hb.1000,
    ncomponents = num_red_dims,
    name = "GLMPCA_approx",
    BSPARAM = BiocSingular::IrlbaParam()
)

spe <- runPCA(
    spe,
    exprs_values = "binomial_deviance_residuals",
    subset_row = hdgs.hb.2000,
    ncomponents = num_red_dims,
    name = "GLMPCA_approx_2000",
    BSPARAM = BiocSingular::IrlbaParam()
)

plotReducedDim(spe, dimred = "GLMPCA_approx", colour_by = "sample_id") 
plotReducedDim(spe, dimred = "GLMPCA_approx_2000", colour_by = "sample_id") 

## Save the processed SPE object

message(Sys.time(), " - Saving HDF5-backed filtered spe")
# spe = saveHDF5SummarizedExperiment(
#     spe, dir = filtered_hdf5_dir, replace = TRUE
# )
# spe = realize(spe)

message(Sys.time(), " - Saving ordinary filtered spe")
# saveRDS(spe, filtered_ordinary_path)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

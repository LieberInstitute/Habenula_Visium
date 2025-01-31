# library("spatialLIBD")
#library("SpatialExperiment")
library("scry")
library("scater")
library("ggplot2"); theme_set(theme_bw())
library("scran")
library("here")
library("tidyverse")

library("BiocParallel")
library("BiocSingular")
library("bluster")
library("PCAtools")

here::here()

## load a filtered spe object
dir_rdata <- here("processed-data", "04_harmony_BayesSpace")
spe <- readRDS(file.path(dir_rdata, "spe_qcED_spatialLIBD_log.rds"))
spe

## Check number of spots in tissue
in_tissue_spots <- sum(as.numeric(map(unique(spe$sample_id), ~ sum(spe$in_tissue[spe$sample_id == .x]))))
print(paste0(" Spots in tissue: ", in_tissue_spots))
map(unique(spe$sample_id), ~ sum(spe$in_tissue[spe$sample_id == .x]))



################################################################################
#   Compute GLM-PCA 
################################################################################

message(Sys.time(), " - Running devianceFeatureSelection()")

# First we computes a deviance statistic for each gene for count data based on a multinomial null model that assumes each feature has a constant rate. 
# Features with large deviance are likely to be informative. Uninformative, low deviance features can be discarded to speed up downstream analyses and reduce memory footprint.

# First we will rank genes based on deviance, to help identify the most biologically informative genes. The actual deviance values are stored in the rowData of the sce/spe object

spe <- devianceFeatureSelection(
  spe,                         
  assay = "counts",  
  fam = "binomial",  # closest approximation to multinomial / Poisson is faster
  nkeep = NULL,      # informative features to be retained
  sorted = FALSE,    # rows sorted in decreasing order of deviance
  batch = as.factor(spe$sample_id)) # batch membership of observations

## Result is a numeric vector containing the deviance statistics for each row.
names(rowData(spe))

spe <- devianceFeatureSelection(spe, assay = "counts", 
                                fam = "poisson", batch = as.factor(spe$sample_id)) 

## A few stats
names(rowData(spe))
binomial_dev <- rowData(spe)$binomial_deviance
poison_dev <- rowData(spe)$poisson_deviance
summary(binomial_dev)
summary(poison_dev)


## plot binomial and poison deviance in first 100 selected genes 

pdf(file.path("plots", "binomial_poisson_dev_visium.pdf"))
par(mfrow = c(2,1))
p1 <- plot(
  sort(rowData(spe)$binomial_deviance, decreasing = TRUE)[1:5000],
  type = "l",
  xlab = "ranked genes",
  ylab = "binomial dev",
  main = "Feature selection with binomial and poisson deviance"
) + abline(v = 1000, lty = 2, col = "red") 
p2 <- plot(
  sort(rowData(spe)$poisson_deviance, decreasing = TRUE)[1:5000],
  type = "l",
  xlab = "ranked genes",
  ylab = "poisson dev",
) + abline(v = 1000, lty = 2, col = "red") 
plts <- p1 / p2
plts
dev.off()

# We can see that the deviance drops sharply after about `x` genes. The remaining genes are probably not informative so we discard them to speed up downstream analysis.

spe2 <- spe[1:1000, ]


# GLM-PCA can reduce the dimensionality of UMI counts to facilitate visualization and/or clustering without needing any normalization.

set.seed(08082024)

## GLM-PCA can reduce the dimensionality of UMI counts to facilitate visualization and/or clustering without needing any normalization.

spe2 <- GLMPCA(spe2, 
               L = 2, 
               assay = "counts",
               minibatch = "stochastic")

# Option 'stochastic' computes a noisy estimate of the full gradient using a random sample of observations at each iteration. 
# Option 'memoized' computes the full data gradient under memory constraints by caching summary  statistics across batches of observations.

reducedDimNames(spe2)
colnames(reducedDim(spe2))


## Visualize the structure

fit <- metadata(spe2)$glmpca
pd <- cbind(as.data.frame(colData(spe2)), fit$factors)
ggplot(pd, aes(x=dim1, y=dim2, colour=sample_id)) + geom_point(size=.8) +
  ggtitle("GLM-PCA applied to high deviance genes")

# colnames(colData(spe2))


# GLM-PCA can be slow for large datasets. A fast approximation is to fit a null model of constant expression for each gene across cells, then fit standard PCA to either the Pearson or deviance residuals from the null model.

## calculate residuals from binomial model

message(Sys.time(), " - Running nullResiduals()")

spe <- nullResiduals( 
  spe,
  assay = "counts",
  fam = "binomial",
  type = "deviance" #,
  # batch = as.factor(spe$sample_id)
)

assayNames(spe)

spe <- nullResiduals(spe, assay="counts", fam ="pearson", type = "deviance")

##  A fast approximation is to fit a null model of constant expression for each gene across cells, then fit standard PCA to either the Pearson or deviance residuals from the null model.

## Get HVDG
hdgs.hb.1000 <-
  rownames(spe)[order(rowData(spe)$binomial_deviance, decreasing = TRUE)][1:1000]


message(Sys.time(), " - Running GLM-PCA")
spe <- runPCA(
  spe,
  exprs_values = "binomial_deviance_residuals",
  subset_row = hdgs.hb.1000,
  ncomponents = 50,
  name = "GLMPCA_approx",
  BSPARAM = BiocSingular::IrlbaParam()
)

plotReducedDim(spe, dimred = "GLMPCA_approx", colour_by = "sample_id")

## The null residuals approach still captures most of the biological structure, but the resolution between clusters is diminished and there is more noise.


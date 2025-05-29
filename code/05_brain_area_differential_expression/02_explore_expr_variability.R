library("here")
library("sessioninfo")
library("SingleCellExperiment")
library("scater")
library("Polychrome")
# copied from https://github.com/LieberInstitute/spatialDLPFC/blob/bd93c980d7653579f81ff1c91c309cea0c7474a6/code/analysis/07_layer_differential_expression/02_explore_expr_variability.R

here::here()

# k <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
args = commandArgs(trailingOnly = TRUE)
k <- as.integer(args[2])

## For testing
## k = 9
if (is.na(k)) {
  k <- 2
}

colors_bayesSpace <- Polychrome::palette36.colors(28)
names(colors_bayesSpace) <- c(1:28)
names(colors_bayesSpace) <-
  paste0("Sp", sprintf("%02d", k), "D", sprintf("%02d", as.integer(names(colors_bayesSpace))))

## output directory
dir_rdata <- here("processed-data", "05_brain_area_differential_expression")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dir_rdata)) ## Check that it was created successfully

dir_plots <- here("plots", "05_brain_area_differential_expression", "02_explore_expr_variability")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dir_plots))

## load spe_pseudo data

spe_pseudo <-
  readRDS(
    file.path(
      dir_rdata,
      paste0("sce_pseudo_PCA_brain_area_k", sprintf("%02d", k), ".rds")
    )
  )

## verify data
dim(spe_pseudo)
dim(reducedDim(spe_pseudo, "PCA"))
rownames(colData(spe_pseudo))
table(spe_pseudo$BayesSpace)
table(spe_pseudo$brain_area2)


# Plot PCs log-transformation; retrieve the PC results

pca_results <- reducedDim(spe_pseudo, "PCA")
head(pca_results)

## Set donor for visualization purposes
spe_pseudo$donor <- spe_pseudo$brain_id

pdf(file = file.path(dir_plots, paste0("sce_pseudo_PC1_k", sprintf("%02d", k), "_main.pdf")), width = 5, height = 5)
plotPCA(spe_pseudo, colour_by = "sample_id", size_by = "sum_umi", shape_by = "donor")
dev.off()

pdf(file = file.path(dir_plots, paste0("sce_pseudo_PC1_k", sprintf("%02d", k), "_sample.pdf")), width = 5, height = 5)
plotPCA(spe_pseudo, colour_by = "sample_id")
dev.off()

pdf(file = file.path(dir_plots, paste0("sce_pseudo_PC1_k", sprintf("%02d", k), "_donor.pdf")), width = 5, height = 5)
plotPCA(spe_pseudo, colour_by = "donor") # , shape_by = "donor"
dev.off()


## Adapted from https://github.com/LieberInstitute/Visium_SPG_AD/blob/6ef1a1225d3dcd115f6272711ab684d050711378/code/11_grey_matter_only/01_create_pseudobulk_data.R#L138-L156

## Explore the resulting data
## - spe it is not just log-transformed, but also normalized by library size (e.g. CPM normalisation)

as.data.frame(colData(spe_pseudo))
pca <- prcomp(t(assays(spe_pseudo)$logcounts))
message(Sys.time(), " % of variance explained for the top 20 PCs:")
metadata(spe_pseudo)
metadata(spe_pseudo) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca)[seq_len(20)])
metadata(spe_pseudo)
pca_pseudo <- pca$x[, seq_len(20)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)

## Already precomputed
# set.seed(20250304)
# spe_pseudo <- scater::runMDS(spe_pseudo, ncomponents = 20)
# spe_pseudo <- scater::runPCA(spe_pseudo, name = "runPCA")


# Calculates the percent of variance explained for first 12 principal components

## Define variables to use
vars <- c(
  "sample_id",
  "donor", 
  "BayesSpace",
  # "age", # it's equivalent to 'donor'
  "brain_area2",
  #"brain_area",
  #"expr_chrM",
  #"expr_chrM_ratio",
  #"nspots",
  #"pmi",
  #"rin",
  "sex"
  #"sum_umi",
)

## Plot PCs with different colors
## Each point here is a sample

# rename levels on colors_bayesSpace variable to match with renamed BayesSpace names

levels(spe_pseudo$BayesSpace)
names(colors_bayesSpace)[1:k]  <- levels(spe_pseudo$BayesSpace)
levels(colors_bayesSpace) <- levels(spe_pseudo$BayesSpace)
levels(colors_bayesSpace)

pdf(file = file.path(dir_plots, paste0("sce_pseudo_PCs_k", sprintf("%02d", k), ".pdf")), width = 8, height = 8)
fontSize = 6

for (var in vars) {
  # var = "BayesSpace"
  ## control legend font.size on BayesSpace plots
  legend.text.font.size <- if (var == "BayesSpace") { if (k < 15) 7 else 5 } else { 7 }
  p <- plotPCA(
    spe_pseudo,
    colour_by = var,
    ncomponents = min(12, length(metadata(spe_pseudo)$PCA_var_explained)),
    point_size = 0.3,
    label_format = c("%s %02i", " (%i%%)"),
    percentVar = metadata(spe_pseudo)$PCA_var_explained
  ) 
  p <- p + ggtitle(paste0("PCA of pseudobulk data with BS k=", as.character(k))) + 
    guides(colour = guide_legend(override.aes = list(size = 6))) +
    theme(text = element_text(size = fontSize), 
          axis.text.y = element_text(size = fontSize),
          axis.text.x = element_text(angle = 90, size = fontSize),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = legend.text.font.size))
  if (var == "BayesSpace") {
    p <- p + scale_color_manual("BayesSpace", values = colors_bayesSpace[1:k])
  }
  print(p)
}
dev.off()

message("PCs plots with different variables done!")

message("Getting variance explained ...")

## Obtain percent of variance explained at the gene level

variance_expl <- scater::getVarianceExplained(spe_pseudo,
                             variables = vars
) 
## Quick inspection
# head(variance_expl)
# summary(variance_expl)


## Now visualize the percent of variance explained across all genes

pdf(file = file.path(dir_plots, paste0("sce_pseudo_gene_explanatory_vars_k", sprintf("%02d", k), ".pdf")))
plotExplanatoryVariables(variance_expl) + ggtitle(paste0("PCA of pseudobulk data with BS k=", as.character(k))) 
dev.off()


message("Plot percent of variance explained across all genes done!")

message("Process completed!!!")



## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


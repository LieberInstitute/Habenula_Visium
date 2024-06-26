## copied from https://github.com/LieberInstitute/Visium_SPG_AD/blob/6ef1a1225d3dcd115f6272711ab684d050711378/code/11_grey_matter_only/01_create_pseudobulk_data.R

# library(sgejobs)
# sgejobs::job_single(
#     name = "01_create_pseudobulk_data",
#     create_shell = TRUE,
#     queue = "shared",
#     memory = "15G",
#     task_num = 28,
#     tc = 10
# )
# To execute the script builder, use: qsub 01_create_pseudobulk_data.sh

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

## For testing
if (FALSE) {
  k <- 2
}

library("here")
library("spatialLIBD")
library("sessioninfo")
library("scater")

## output directory
# dir_rdata <- here::here(
#   "processed-data",
#   "rdata",
#   "spe",
#   "07_layer_differential_expression"
# )
dir_rdata <- here("processed-data", "05_layer_differential_expression")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dir_rdata)) ## Check that it was created successfully

## load spe data
# load(
#   here(
#     "processed-data",
#     "rdata",
#     "spe",
#     "01_build_spe",
#     "spe_filtered_final_with_clusters.Rdata"
#   ),
#   verbose = TRUE
# )
spe_in <- here("processed-data", "04_harmony_BayesSpace", "spe_harmony.rds")
spe <- readRDS(spe_in)

## Import BayesSpace clusters
colnames(colData(spe))
clusters_BayesSpace_dir <- here("processed-data", "04_harmony_BayesSpace", "clusters_BayesSpace")
spe <- cluster_import(spe,
                      cluster_dir = clusters_BayesSpace_dir,              
                      prefix = ""
)

## Convert from character to a factor

# Quick inspection
colData(spe)[grep("BayesSpace_harmony", colnames(colData(spe)))]
#length(grep("BayesSpace_harmony", colnames(colData(spe))))

spe$BayesSpace <-
  factor(
    paste0("Sp", sprintf("%02d", k), "D", sprintf("%02d", colData(spe)[[paste0("BayesSpace_harmony_k", sprintf("%02d", k))]]))
  )

## pseudobulk across a given BayesSpace k
sce_pseudo <-
  registration_pseudobulk(spe,
                          var_registration = "BayesSpace",
                          var_sample_id = "sample_id",
                          min_ncells = 10
  )
dim(sce_pseudo)
# [1] 13508     15 with 5 capture areas and k=5

## Rename "region" into "position" for consistency with
## https://github.com/LieberInstitute/DLPFC_snRNAseq
#sce_pseudo$position <- sce_pseudo$region

## Simplify the colData()  for the pseudo-bulked data
colnames(colData(sce_pseudo))
colData(sce_pseudo) <- colData(sce_pseudo)[, sort(c(
  "age",
  "sample_id",
  "BayesSpace",
  "subject",
  "sex",
  # "position",
  "diagnosis",
  "ncells"
))]

## Explore the resulting data
options(width = 400)
as.data.frame(colData(sce_pseudo))

## Compute PCs
## Adapted from https://github.com/LieberInstitute/spatialDLPFC/blob/f47daafa19b02e6208c7e0a9bc068367f806206c/code/analysis/09_region_differential_expression/preliminary_analysis.R#L60-L68

pca <- prcomp(t(assays(sce_pseudo)$logcounts))
# length(pca$sdev) #number of components processed 
# summary(pca) 
# Importance of components:
#                         PC1     PC2     PC3      PC4      PC5      PC6      PC7      PC8      PC9    PC10    PC11    PC12   PC13    PC14      PC15
# Standard deviation     60.90 36.1595 29.9354 20.08013 19.07484 16.91213 12.60943 12.20756 11.23785 9.15730 8.26029 7.97064 7.7150 7.17131 5.309e-14
# Proportion of Variance  0.48  0.1692  0.1160  0.05218  0.04709  0.03702  0.02058  0.01929  0.01634 0.01085 0.00883 0.00822 0.0077 0.00666 0.000e+00
# Cumulative Proportion   0.48  0.6493  0.7652  0.81742  0.86451  0.90153  0.92211  0.94139  0.95774 0.96859 0.97742 0.98564 0.9933 1.00000 1.000e+00
# print(pca[1])
# plot(pca, "PCA of pseudobulk data")
# plot(pca$x[,1],pca$x[,2])
# biplot(pca)

# Set number of components equal to pseudo bulk groups. Avoid an error triggered when <20 pseudo bulk groups are processed by default.
if ((n_components <- length(pca$sdev)) > 21) { n_components <- 20 } 

message(Sys.time(), " % of variance explained for the top ", n_components ," PCs:")
metadata(sce_pseudo)
metadata(sce_pseudo) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca))[seq_len(n_components)] #20
pca_pseudo <- pca$x[, seq_len(n_components)] #20
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(sce_pseudo) <- list(PCA = pca_pseudo)

## Compute some reduced dims
set.seed(20240626)
# sce_pseudo <- scater::runMDS(sce_pseudo, ncomponents = (n_components-1)) #20
# sce_pseudo <- scater::runPCA(sce_pseudo, name = "runPCA")

## We don't want to model the pathology groups as integers / numeric
## so let's double check this
stopifnot(is.factor(sce_pseudo$BayesSpace))

## For the spatialLIBD shiny app
rowData(sce_pseudo)$gene_search <-
  paste0(
    rowData(sce_pseudo)$gene_name,
    "; ",
    rowData(sce_pseudo)$gene_id
  )

## Load pathology colors
## This info is used by spatialLIBD v1.7.18 or newer
# source(here("code", "analysis", "colors_bayesSpace.R"), echo = TRUE, max.deparse.length = 500)
# names(colors_bayesSpace) <-
#   paste0("Sp", sprintf("%02d", k), "D", sprintf("%02d", as.integer(names(colors_bayesSpace))))
# sce_pseudo$BayesSpace_colors <- colors_bayesSpace[as.character(sce_pseudo$BayesSpace)]

## save RDS file
saveRDS(
  sce_pseudo,
  file = file.path(
    dir_rdata,
    paste0("sce_pseudo_BayesSpace_k", sprintf("%02d", k), ".rds")
  )
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
## copied from https://github.com/LieberInstitute/Visium_SPG_AD/blob/6ef1a1225d3dcd115f6272711ab684d050711378/code/11_grey_matter_only/01_create_pseudobulk_data.R

# library(slurmjobs)
# slurmjobs::job_single('01_create_pseudobulk_data', create_shell = TRUE, memory = '20G', command = "01_create_pseudobulk_data.R")
# To submit the job use: sbatch 01_create_pseudobulk_data.sh

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

spe$BayesSpace <- factor(
    paste0("Sp", sprintf("%02d", k), "D", sprintf("%02d", colData(spe)[[paste0("BayesSpace_harmony_k", sprintf("%02d", k))]]))
)
# head(unique(spe$BayesSpace))

## pseudobulk across a given BayesSpace k
sce_pseudo <-
  registration_pseudobulk(spe,
                          var_registration = "BayesSpace",
                          var_sample_id = "sample_id",
                          min_ncells = 10
  )
dim(sce_pseudo)

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

## Explore pca
# length(pca$sdev) 
# summary(pca) 
# print(pca[1])
# plot(pca, paste0("PCA of pseudobulk data with BS k=", as.character(k)))
# plot(pca$x[,1],pca$x[,2])
# biplot(pca)


# Set number of components equal to pseudo bulk groups. Avoid an error triggered when number of components <20 pseudo-bulked groups. Default componenets = 20
if ((n_components <- length(pca$sdev)) > 21) { n_components <- 20 } 

message(Sys.time(), " % of variance explained for the top ", n_components ," PCs:")
metadata(sce_pseudo)
# metadata(sce_pseudo) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca))[seq_len(n_components)]
metadata(sce_pseudo) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca)[seq_len(n_components)])
pca_pseudo <- pca$x[, seq_len(n_components)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(sce_pseudo) <- list(PCA = pca_pseudo)
# plotPCA(sce_pseudo, colour_by = "sample_id", ncomponents = n_components, point_size = 1) 

## Compute some reduced dims
set.seed(20240626)
sce_pseudo <- scater::runMDS(sce_pseudo, ncomponents = (n_components-1)) #20
sce_pseudo <- scater::runPCA(sce_pseudo, name = "runPCA")
# Warning in (function (A, nv = 5, nu = nv, maxit = 1000, work = nv + 7, reorth = TRUE,  :
#                         You're computing too large a percentage of total singular values, use a standard svd instead.

## Double check the BayesSpace meta are factors
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
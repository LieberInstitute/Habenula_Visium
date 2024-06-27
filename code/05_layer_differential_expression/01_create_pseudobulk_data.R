## copied from https://github.com/LieberInstitute/Visium_SPG_AD/blob/6ef1a1225d3dcd115f6272711ab684d050711378/code/11_grey_matter_only/01_create_pseudobulk_data.R

# library(slurmjobs)
# slurmjobs::job_single('01_create_pseudobulk_data', create_shell = TRUE, memory = '20G', command = "01_create_pseudobulk_data.R")
# To submit the job use: sbatch 01_create_pseudobulk_data.sh

k <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

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
head(spe$key[1:5])
spe <- cluster_import(spe,
                      cluster_dir = clusters_BayesSpace_dir,              
                      prefix = "",
                      overwrite = TRUE
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

message('Pseudobulk completed ')


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
  "diagnosis",
  "ncells"
))]

## Explore the resulting data
options(width = 400)
as.data.frame(colData(sce_pseudo))

## Compute PCs
## Adapted from https://github.com/LieberInstitute/spatialDLPFC/blob/f47daafa19b02e6208c7e0a9bc068367f806206c/code/analysis/09_region_differential_expression/preliminary_analysis.R#L60-L68

message('Processing PCA')

pca <- prcomp(t(assays(sce_pseudo)$logcounts))
dim(pca$x)
## Explore pca
# length(pca$sdev) 
# summary(pca) 
# print(pca$x)
# plot(pca, paste0("PCA of pseudobulk data with BS k=", as.character(k)))
# plot(pca$x[,1],pca$x[,2])

message(Sys.time(), " % of variance explained for the top PCs:")
metadata(sce_pseudo)
metadata(sce_pseudo) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca)) #[seq_len(20)])
metadata(sce_pseudo)
#pca_pseudo <- pca$x[, seq_len(n_components)]
colnames(pca$x) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca$x))))
head(pca$x)
reducedDims(sce_pseudo) <- list(PCA = pca$x)
#plotPCA(sce_pseudo, colour_by = "sample_id", 2, point_size = 1) 

message(' Compute reduced dims')

## Compute some reduced dims
set.seed(20240626)
#sce_pseudo <- scater::runMDS(sce_pseudo, ncomponents = 2) #20
#sce_pseudo <- scater::runPCA(sce_pseudo, name = "runPCA")

## Double check the BayesSpace meta are factors
stopifnot(is.factor(sce_pseudo$BayesSpace))

# ## For the spatialLIBD shiny app
# rowData(sce_pseudo)$gene_search <-
#   paste0(
#     rowData(sce_pseudo)$gene_name,
#     "; ",
#     rowData(sce_pseudo)$gene_id
#   )

message(' Saving pseudobulk with reduced dims for BS k=', k)

## save RDS file
saveRDS(
  sce_pseudo,
  file = file.path(
    dir_rdata,
    paste0("sce_pseudo_BayesSpace_k", sprintf("%02d", k), ".rds")
  )
)

message(' Process completed!')

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

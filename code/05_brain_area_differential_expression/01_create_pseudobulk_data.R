## copied from https://github.com/LieberInstitute/Visium_SPG_AD/blob/6ef1a1225d3dcd115f6272711ab684d050711378/code/11_grey_matter_only/01_create_pseudobulk_data.R

# library(slurmjobs)
# slurmjobs::job_single('01_create_pseudobulk_data', 
#                       create_shell = TRUE, memory = '60G', 
#                       command = "01_create_pseudobulk_data.R", 
#                       partion = "katun")


k <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## For testing
if (is.na(k)) {
  k <- 2
}

library("here")
library("spatialLIBD")
library("ggplot2")
library("gridExtra")
library("sessioninfo")
library("scater")

dir_rdata <- here("processed-data", "05_brain_area_differential_expression")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dir_rdata)) ## Check that it was created successfully

## load spe data
spe_in <- here("processed-data", "04_harmony_BayesSpace", "spe_harmony.rds")
spe <- readRDS(spe_in)

## Import BayesSpace clusters
clusters_BayesSpace_dir <- here("processed-data", "04_harmony_BayesSpace", "clusters_BayesSpace")
# head(spe$key[1:5])
spe <- cluster_import(spe,
                      cluster_dir = clusters_BayesSpace_dir,              
                      prefix = "",
                      overwrite = TRUE
)
# Overwriting 'spe$key'. Set 'overwrite = FALSE' if you do not want to overwrite it.

## Convert from character to a factor

# Quick inspection
colData(spe)[grep("BayesSpace_harmony", colnames(colData(spe)))]
#length(grep("BayesSpace_harmony", colnames(colData(spe))))

spe$BayesSpace <- factor(
    paste0("Sp", sprintf("%02d", k), "D", sprintf("%02d", colData(spe)[[paste0("BayesSpace_harmony_k", sprintf("%02d", k))]]))
)
# head(unique(spe$BayesSpace))
# [1] Sp02D02 Sp02D01
# Levels: Sp02D01 Sp02D02

# Add a new column based on brain_id condition - This will be used as variable for registration

#colnames(colData(spe))
table(colData(spe)$brain_area)
colData(spe)$brain_area_DEG <- ifelse(colData(spe)$brain_area == "AR6" | colData(spe)$brain_area == "AL5", "Anterior", "Posterior")
table(colData(spe)$brain_area_DEG)
# Anterior Posterior 
# 22571     10838
# Convert relevant variables
table(spe_pseudo$brain_area_DEG)
# Anterior Posterior
#   7         5
spe_pseudo$brain_area_DEG <- factor(spe_pseudo$brain_area_DEG, levels = c("Anterior", "Posterior"))
levels(spe_pseudo$brain_area_DEG)
table(colData(spe)$brain_id)

spe_pseudo <-
  registration_pseudobulk(
    spe,
    var_registration = "brain_area_DEG",
    var_sample_id = "sample_id",
    covars = "brain_id",
    min_ncells = 10
  )

## Drop unused var_registration levels if we had to drop some due to min_nspots:

## drop levels not used
spe_pseudo$brain_area_DEG <- droplevels(spe_pseudo$brain_area_DEG)
spe_pseudo$brain_area_DEG <- factor(spe_pseudo$brain_area_DEG, levels = c("Anterior", "Posterior"))
levels(spe_pseudo$brain_area_DEG)
table(spe_pseudo$brain_area_DEG)
## set numeric to avoid error reading age variable
spe_pseudo$age <- as.numeric(spe_pseudo$age)

message('Levels unused on pseudobulk `brain_area_DEG` dropped ')

message('Pseudobulk completed ')

# # Compute mitochondrial expression ratio
# is_mito <- which(seqnames(spe_pseudo) == "chrM")
# spe_pseudo$expr_chrM <- colSums(counts(spe_pseudo)[is_mito, , drop = FALSE])
# spe_pseudo$sum_umi <- colSums(counts(spe_pseudo))
# spe_pseudo$expr_chrM_ratio <- spe_pseudo$expr_chrM / spe_pseudo$sum_umi


## Simplify the colData()  for the pseudo-bulked data
colnames(colData(spe_pseudo))
colData(spe_pseudo) <- colData(spe_pseudo)[, sort(c(
  "age",
  "sample_id",
  # "BayesSpace",
  "brain_area_DEG",
  "brain_id", # equivalent to subject / donor / ethnicity
  "sex",
  "diagnosis",
  "ncells"
))]

## Explore the resulting data
options(width = 400)
as.data.frame(colData(spe_pseudo))

## Compute PCs
## Adapted from https://github.com/LieberInstitute/spatialDLPFC/blob/f47daafa19b02e6208c7e0a9bc068367f806206c/code/analysis/09_region_differential_expression/preliminary_analysis.R#L60-L68

message('Processing PCA')

# First, performed PCA manually using prcomp()

max_components <- min(dim(spe_pseudo)) - 1
print(max_components)
n_components <- min(n_components, max_components)
pca <- prcomp(t(assays(spe_pseudo)$logcounts), center = TRUE, scale. = TRUE)
dim(pca$x)
names(pca)
# Store PCA coordinates
reducedDims(spe_pseudo)$PCA <- pca_result$x
reducedDims(spe_pseudo)

# Set number of components equal to pseudo bulk groups. Avoid an error triggered when n_components <20 pseudo bulk groups. Default=20.
n_components <- length(pca$sdev)
if (n_components > 21) { n_components <- 20 }

message(Sys.time(), " % of variance explained for the top ", n_components ," PCs:")
metadata(spe_pseudo) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca)[seq_len(n_components)]) #[seq_len(20)])
# metadata(spe_pseudo)
colnames(pca$x) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca$x))))
# head(pca$x)
# View PCA values
head(reducedDim(spe_pseudo, "PCA"))
reducedDims(spe_pseudo) <- list(PCA = pca$x)
colnames(colData(spe_pseudo))

# Quick inspection 
plt1 = plotPCA(spe_pseudo, colour_by = "sample_id", n_components, point_size = 3) 
plt2 =  plotPCA(spe_pseudo, colour_by = "brain_area_DEG", n_components, point_size = 3) 
plt = grid.arrange(plt1, plt2, ncol=2)

## Compute some reduced dims
message('/nProcessing MDS and scarter runPCA')

set.seed(20240626)
spe_pseudo <- scater::runMDS(spe_pseudo, name = "runMDS", ncomponents = (n_components-1))
# spe_pseudo <- scater::runPCA(spe_pseudo, name = "runPCA", ncomponents = n_components) 

## Double check the brain_area_DEG meta are factors
stopifnot(is.factor(spe_pseudo$brain_area_DEG))

## For the spatialLIBD shiny app
rowData(spe_pseudo)$gene_search <-
  paste0(
    rowData(spe_pseudo)$gene_name,
    "; ",
    rowData(spe_pseudo)$gene_id
  )

message(' Saving pseudobulk with reduced dims for BS k=', k, " using brain_area as var for registration")

## save RDS file
saveRDS(
  spe_pseudo,
  file = file.path(
    dir_rdata,
    paste0("sce_pseudo_brain_area_", sprintf("%02d", k), ".rds")
  )
)


message(' Process completed!')

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

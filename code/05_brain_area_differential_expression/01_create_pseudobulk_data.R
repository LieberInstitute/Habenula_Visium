library("here")
library("spatialLIBD")
library("tidyverse")
library("ggplot2")
library("gridExtra")
library("sessioninfo")
library("scater")
library("BiocSingular") # Force svd method on pca
library("compositions")
#install.packages("compositions")

## copied from https://github.com/LieberInstitute/Visium_SPG_AD/blob/6ef1a1225d3dcd115f6272711ab684d050711378/code/11_grey_matter_only/01_create_pseudobulk_data.R

# k <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
args = commandArgs(trailingOnly = TRUE)
k <- as.integer(args[2])

## For testing
if (is.na(k)) {
  k <- 2
}

message("Processing pseudobullk for k = ", k)

dir_rdata <- here("processed-data", "05_brain_area_differential_expression")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dir_rdata)) ## Check that it was created successfully
dir_csv <- here(
  "processed-data",
  "05_brain_area_differential_expression",
  "stats_summary_csv"
)
dir.create(dir_csv, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dir_csv))

## load spe data

spe_in <- here("processed-data", "04_harmony_BayesSpace", "spe_harmony.rds")
spe <- readRDS(spe_in)

## Import BayesSpace clusters

clusters_BayesSpace_dir <- here(
  "processed-data",
  "04_harmony_BayesSpace",
  "clusters_BayesSpace"
)
# head(spe$key[1:5])
spe <- cluster_import(
  spe,
  cluster_dir = clusters_BayesSpace_dir,
  prefix = "",
  overwrite = TRUE
)
# Overwriting 'spe$key'. Set 'overwrite = FALSE' if you do not want to overwrite it.

## Prepare data to pseudobulk

# Quick inspection
#colData(spe)[grep("BayesSpace_harmony", colnames(colData(spe)))]
#length(grep("BayesSpace_harmony", colnames(colData(spe))))

k_nice <- sprintf("%02d", k)

spe$BayesSpace <- factor(
  paste0(
    "Sp",
    k_nice,
    "D",
    sprintf(
      "%02d",
      colData(spe)[[paste0("BayesSpace_harmony_k", sprintf("%02d", k))]]
    )
  )
)
head(unique(spe$BayesSpace))
# [1] Sp02D02 Sp02D01
# Levels: Sp02D01 Sp02D02

levels(colData(spe)$BayesSpace)
# [1] "Sp02D01" "Sp02D02" 

## load spatial-registration to pull short annotated labels

dir_labels <- here("processed-data", 
                  "06_spatial_registration_vs_snRNA-seq", 
                  "cor_BayesSpace_vs_snRNA-seq_top100_Hb_merged.Rdata")
load(dir_labels)
# [3] "cor_broad"               "cor_fine"
levels(cor_broad)
# sort(rownames(cor_broad[[k-1]]))

# rename the levels to make id readable in the plots

colData(spe)$BayesSpace <- factor(colData(spe)$BayesSpace)
levels(colData(spe)$BayesSpace) <- c(sort(rownames(cor_broad[[k-1]])))
levels(colData(spe)$BayesSpace)
# [1] "Sp02D01 ~ Oligo"      "Sp02D02 ~ Inhib.Thal"

# Add a new column based on brain_id condition - This will be used as variable for exploring variation

#colnames(colData(spe))
table(colData(spe)$sample_id)
table(colData(spe)$brain_area)
table(colData(spe)$BayesSpace)

## Assign new 'brain-area' based in posterior-anterior locations defined by KDM based on RNAScope

colData(spe)$brain_area2 <- case_when(
  colData(spe)$sample_id == "V13B23-285_A1" ~ "G0",
  colData(spe)$sample_id == "V13B23-285_B1" |
    colData(spe)$sample_id == "V14F07-340_A1" |
    colData(spe)$sample_id == "V13B23-280_A1" ~
    "G1",
  colData(spe)$sample_id == "V13B23-285_C1" |
    colData(spe)$sample_id == "V14F07-340_B1" |
    colData(spe)$sample_id == "V13B23-280_B1" ~
    "G2",
  colData(spe)$sample_id == "V13B23-285_D1" |
    colData(spe)$sample_id == "V14F07-340_C1" |
    colData(spe)$sample_id == "V13B23-280_C1" ~
    "G3",
  colData(spe)$sample_id == "V14F07-340_D1" |
    colData(spe)$sample_id == "V13B23-280_D1" ~
    "G4"
)
table(colData(spe)$brain_id, colData(spe)$brain_area2)
#           G0   G1   G2   G3   G4
# Br8518 3803 3119 3201 3118    0
# Br9037    0 3117 2921 3580 3515
# Br9090    0 1669 1689 1794 1883

# Previous definition
# colData(spe)$brain_area_DEG <- ifelse(
#   colData(spe)$brain_area == "AR6" | colData(spe)$brain_area == "AL5",
#   "Anterior",
#   "Posterior"
# )
# table(colData(spe)$brain_area_DEG)
# Anterior Posterior
# 22571     10838

table(colData(spe)$brain_id)
# Br8518 Br9037 Br9090
# 13241  13133   7035

############################

message("Processing BayesSpace k=", k_nice)

# Perform pseudobulk across BayesSpace and sample_id

spe_pseudo_k <- aggregateAcrossCells(
  spe,
  DataFrame(
    BayesSpace = spe[[paste0("BayesSpace_harmony_k", k_nice)]],
    reg_sample_id = spe$sample_id
  )
)
levels(spe_pseudo_k$BayesSpace)

message("Aggregation completed for k=", k_nice)
message(
  "Dimensions of summed data: ",
  paste(dim(spe_pseudo_k), collapse = " x ")
)

# Rename ncells to nspots

colData(spe_pseudo_k)$nspots <- colData(spe_pseudo_k)$ncells
colData(spe_pseudo_k)$ncells <- NULL # Remove the old column
colnames(spe_pseudo_k) <- spe_pseudo_k$sample_id


# Exploring nspots

min_nspots <- 10
message("Total nspots: ", sum(spe_pseudo_k$nspots))
message(
  "Number of groups with nspots < ",
  min_nspots,
  ": ",
  sum(spe_pseudo_k$nspots < min_nspots)
)
message("Summary of nspots:")
print(summary(spe_pseudo_k$nspots))

## Adapted from https://github.com/LieberInstitute/spatialLIBD/blob/devel/R/registration_pseudobulk.R#L137-L154
## Drop pseudo-bulked samples that had low initial contribution of raw-samples.
## That is, pseudo-bulked samples that are not benefiting from the pseudo-bulking process to obtain higher counts.
if (!is.null(min_nspots)) {
  message(
    Sys.time(),
    " dropping ",
    sum(spe_pseudo_k$nspots < min_nspots),
    " pseudo-bulked samples that are below 'min_nspots'."
  )
  spe_pseudo_k <- spe_pseudo_k[, spe_pseudo_k$nspots >= min_nspots]
}

# Compute mitochondrial expression ratio

is_mito <- which(seqnames(spe_pseudo_k) == "chrM")
spe_pseudo_k$expr_chrM <- colSums(counts(spe_pseudo_k)[is_mito, , drop = FALSE])
spe_pseudo_k$sum_umi <- colSums(counts(spe_pseudo_k))
spe_pseudo_k$expr_chrM_ratio <- spe_pseudo_k$expr_chrM / spe_pseudo_k$sum_umi

spe_pseudo_k$age <- as.numeric(spe_pseudo_k$age)

# Convert relevant variables to factors

if (is.factor(spe_pseudo_k$BayesSpace)) {
  ## Drop unused var_registration levels if we had to drop some due to min_nspots:
  ## registration_variable equivalent here: BayesSpace_pseudo
  spe_pseudo_k$BayesSpace <- droplevels(spe_pseudo_k$BayesSpace)
  levels(spe_pseudo_k$BayesSpace)
}

spe_pseudo_k$brain_area2 <- factor(spe_pseudo_k$brain_area2)
if (is.factor(spe_pseudo_k$BayesSpace_pseudo)) {
  spe_pseudo_k$brain_area2 <- droplevels(spe_pseudo_k$brain_area2)
}
levels(spe_pseudo_k$brain_area2)

## Compute the logcounts

message(Sys.time(), " normalize expression")

assays(spe_pseudo_k)

logcounts(spe_pseudo_k) <-
  edgeR::cpm(edgeR::calcNormFactors(spe_pseudo_k), log = TRUE, prior.count = 1)
# rownames(assays(spe_pseudo_k)$logcounts)
# colnames(assays(spe_pseudo_k)$logcounts)
dim(reducedDim(spe_pseudo_k))
# [1] 24 10

# # calculate the number of cells per (sample_id + BayesSpace cluster)
# Adapted from: https://github.com/LieberInstitute/dlpfc_asd/blob/2b83eeb9572bd7d37505e8db6e20bb3ded09c2c1/code/06_differential_expression/01_create_pseudobulk_data.R#L129

message("-------------------------------------------------------------")
# Get BayesSpace cluster assignments for the current k
cluster_ids <- spe[[paste0("BayesSpace_harmony_k", k_nice)]]
# Create a data frame with sample_id and cluster_ids
df <- data.frame(sample_id = spe$sample_id, cluster_ids = cluster_ids)
# Create a table between sample_id and cluster_ids
k_table <- table(df$sample_id, df$cluster_ids)
# Add a 'Total' column to the table by summing across rows (sum of cells for each sample_id)
k_table <- cbind(k_table, Total = rowSums(k_table))
# Add a ncells from spe_pseudo_k
k_table <- cbind(k_table, nspots = spe_pseudo_k$nspots)
k_table_subset <- k_table[
  colnames(spe_pseudo_k),
  grepl("^[[:digit:]]+$", colnames(k_table))
]

## Compute ILR
k_table_ilr <- ilr(k_table_subset)
## Note that this is basically the same as
## ilr(k_table_subset / rowSums(k_table_subset))
##rowSums(k_table_subset / rowSums(k_table_subset))  equal to 1
colnames(k_table_ilr) <- paste0(
  "ILR_SpD",
  k_nice,
  "_",
  seq_len(ncol(k_table_ilr))
)
colData(spe_pseudo_k) <- cbind(
  colData(spe_pseudo_k),
  k_table_subset,
  as.data.frame(k_table_ilr)
)

# Print the table of cell counts per sample_id and cluster
message("Cell counts per sample_id and cluster for k=", k_nice)
print(k_table)

message("-------------------------------------------------------------")

## save table with basic stats
write.csv(
  k_table,
  row.names = TRUE,
  quote = FALSE,
  here(
    dir_csv,
    paste0("k", sprintf("%02d", k), "_basic_stats.csv")
  )
)

message('Pseudobulk completed ')


## Simplify the colData()  for the pseudo-bulked data

# colnames(colData(spe_pseudo_k))
colData(spe_pseudo_k) <- colData(spe_pseudo_k)[, sort(c(
  "sample_id",
  "brain_id",
  "age",
  "sex",
  "diagnosis",
  "brain_area",
  "brain_area2",
  "nspots",
  "sum_umi",
  "expr_chrM",
  "expr_chrM_ratio",
  "rin",
  "pmi",
  "BayesSpace"
))]

## Explore the resulting data

options(width = 400)
sp_table = as.data.frame(colData(spe_pseudo_k))
## save table with general information about spatial-domains
write.csv(
  sp_table,
  row.names = TRUE,
  quote = FALSE,
  here(
    dir_csv,
    paste0("k", sprintf("%02d", k), "_SpatialD_info.csv")
  )
)


message('Processing PCA')

set.seed(01042025)

spe_pseudo <- spe_pseudo_k

## Compute some reduced dims
message('Processing MDS and scarter runPCA')

# keep ncomponents below the number of cells
ncomponents = min(50, ncol(spe_pseudo) - 1, nrow(spe_pseudo) - 1)

spe_pseudo <- scater::runMDS(
  spe_pseudo,
  name = "runMDS",
  ncomponents = ncomponents
)


# Fix error when computing runPCA with irlba  - scater default -  too many components relative to the number of features
# -- force exact SVD passing a BiocSingularParam object

spe_pseudo <- scater::runPCA(
  spe_pseudo,
  ncomponents = ncomponents,
  BSPARAM = ExactParam()
)

dim(reducedDim(spe_pseudo, "PCA"))
# [1] 24 23

# message(
#   Sys.time(),
#   " % of variance explained for the top ",
#   n_components,
#   " PCs:"
# )
# ##  computes the percent of variance explained by each of the principal components - created with prcomp
# metadata(spe_pseudo) <- list(
#   "PCA_var_explained" = jaffelab::getPcaVars(pca)[seq_len(n_components)]
# )
# # metadata(spe_pseudo)
# colnames(pca$x) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca$x))))
# # head(pca$x)
# # View PCA values
# head(reducedDim(spe_pseudo, "PCA"))
# reducedDims(spe_pseudo) <- list(PCA = pca$x)

## For the spatialLIBD shiny app

rowData(spe_pseudo)$gene_search <-
  paste0(
    rowData(spe_pseudo)$gene_name,
    "; ",
    rowData(spe_pseudo)$gene_id
  )

message(
  ' Saving pseudobulk with reduced dims for BS k=',
  k,
  " using brain_area as var for registration"
)

## save RDS file
saveRDS(
  spe_pseudo,
  file = file.path(
    dir_rdata,
    paste0("sce_pseudo_PCA_brain_area_k", sprintf("%02d", k), ".rds")
  )
)


message(' Process completed!')


# library(slurmjobs)
# slurmjobs::job_single('01_create_pseudobulk_data',
#                       create_shell = TRUE, memory = '60G',
#                       command = "01_create_pseudobulk_data.R",
#                       partion = "katun")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

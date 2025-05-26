########################################################################
## copied from https://github.com/LieberInstitute/Visium_SPG_AD/blob/6ef1a1225d3dcd115f6272711ab684d050711378/code/11_grey_matter_only/01_create_pseudobulk_data.R
## Implementation: CSC
## Data: XX
## For 60 to 80k spots: $srun --pty --mem=30GB --x11 bash
########################################################################


library("here")
library("spatialLIBD")
library("tidyverse")
library("ggplot2")
library("gridExtra")
library("sessioninfo")
library("scater")
library("BiocSingular") # Force svd method on pca
library("compositions")


k <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## For testing
if (is.na(k)) {
  k <- 2
}

k_nice <- sprintf("%02d", k)

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
# class: SpatialExperiment 
# dim: 26610 33409 

## Import BayesSpace clusters
colnames(colData(spe))
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


message("Processing pseudobullk for k = ", k_nice)


## Prepare data to pseudobulk

# Quick inspection
#colData(spe)[grep("BayesSpace_harmony", colnames(colData(spe)))]
#length(grep("BayesSpace_harmony", colnames(colData(spe))))
colnames(colData(spe))[grepl("BayesSpace_harmony*", colnames(colData(spe)))]
table(spe$brain_id)
table(spe$sample_id)

## set `BayesSpace` as factor
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

levels(colData(spe)$BayesSpace)
# [1] "Sp02D01" "Sp02D02"

# ## load annotated data at 'Broad' level to extract short annotated labels -- need to be check (CSC)
# 
# dir_labels <- here("processed-data",
#                   "06_spatial_registration_vs_snRNA-seq",
#                   "cor_BayesSpace_vs_snRNA-seq_top100.Rdata")
# load(dir_labels)
# # [3] "cor_broad"               "cor_fine"
# 
# # rename the levels to make them readable in the plots (uses Hb pilot Broad annotations)
# 
# levels(colData(spe)$BayesSpace) <- c(sort(rownames(cor_broad[[k-1]])))
# levels(colData(spe)$BayesSpace)
# # [1] "Sp02D01 ~ Oligo"      "Sp02D02 ~ Inhib.Thal"
# 
# ## we need to double-check if the annotated BayesSpace names/levels do have syntactically valid names
# #  - this avoid error further when computing other process. g.e: 'registration_stats_pairwise'
# levels(spe$BayesSpace) <- gsub("\\s*~\\s*", ".", levels(spe$BayesSpace))  # Replace ~
# levels(spe$BayesSpace) <- gsub("/", "_", levels(spe$BayesSpace))          # Replace /
# levels(spe$BayesSpace) <- gsub("\\*$", "", levels(spe$BayesSpace))        # Remove *
# # Check
# levels(spe$BayesSpace)
# head(spe$BayesSpace)


## quick inspection: check how many genes expressed by cluster we have before pseudobulk

table(colData(spe)$brain_id, colData(spe)$sample_id)
table(colData(spe)$brain_id, colData(spe)$BayesSpace)

# number of genes before pseudobulk
expr_mat <- assay(spe, "logcounts")
clusters <- spe$BayesSpace
unique_clusters <- levels(factor(clusters))
genes_per_cluster <- sapply(unique_clusters, function(clust) {
    # Subset expression matrix to spots in this cluster
    cluster_expr <- expr_mat[, clusters == clust]
    # Count genes with at least one non-zero value in the cluster
    sum(rowSums(cluster_expr > 0) > 0)
})
gene_counts_df <- data.frame(
    Cluster = unique_clusters,
    Num_Expressed_Genes = genes_per_cluster,
    row.names = NULL
)
message("Number of expressed genes before pseudobulk:  ")
gene_counts_df
#               Cluster Num_Expressed_Genes
# 1      Sp03D01.Oligo               18798
# 2   Sp03D02.Habenula               26223
# 3 Sp03D03.Inhib.Thal               21454


## Assign new column 'brain-area2' based on `anterior to posterior` regions defined by KDM based on RNAScope

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

# Make sure it's a factor
spe$brain_area2 <- factor(spe$brain_area2)
table(spe$BayesSpace, spe$brain_area2)
#                       G0   G1   G2   G3   G4
# Sp03D01.Oligo       700 1538 1231 1335 1154
# Sp03D02.Habenula   3100 6031 6239 6763 3963
# Sp03D03.Inhib.Thal    3  336  341  394  281

table(colData(spe)$brain_id, colData(spe)$brain_area2)
#           G0   G1   G2   G3   G4
# Br8518 3803 3119 3201 3118    0
# Br9037    0 3117 2921 3580 3515
# Br9090    0 1669 1689 1794 1883

table(colData(spe)$brain_id)
# Br8518 Br9037 Br9090
# 13241  13133   7035


## save new spe object with BayesSpace clusters

spe_in <- here("processed-data", "05_brain_area_differential_expression", "spe_harmony_ann.rds")
saveRDS(
  spe, file = file.path(spe_in)
)


############################

message("Processing pseudobulk for BayesSpace k=", k_nice)

spe_pseudo_k <-
    registration_pseudobulk(
        spe,
        var_registration = "BayesSpace",
        var_sample_id = "sample_id",
        covars = "brain_id",
        min_ncells = 10
    )
dim(spe_pseudo_k)

#colnames(colData(spe_pseudo_k))
## list domains created
rownames(colData(spe_pseudo_k))
table(spe_pseudo_k$sample_id)
table(spe_pseudo_k$brain_id)

## drop levels not used
table(spe_pseudo_k$BayesSpace)
unique(spe_pseudo_k$BayesSpace)
levels(spe_pseudo_k$BayesSpace)
spe_pseudo_k$BayesSpace <- droplevels(spe_pseudo_k$BayesSpace)
levels(spe_pseudo_k$BayesSpace)
table(spe_pseudo_k$BayesSpace)

message('Levels unused on pseudobulk `BayesSpace` dropped ')

message('Pseudobulk completed ')

## quick inspection: check how many genes expressed by cluster we have after pseudobulk

expr_mat <- assay(spe_pseudo_k, "counts")  # or "counts" if logcounts not available
clusters <- spe_pseudo_k$BayesSpace
unique_clusters <- unique(clusters)
# Count expressed genes per cluster
genes_per_cluster <- sapply(unique_clusters, function(clust) {
    cluster_expr <- expr_mat[, clusters == clust]
    # Count genes with expression > 0 in at least one pseudobulked sample
    sum(rowSums(cluster_expr > 0) > 0)
})
gene_counts_df <- data.frame(
    Cluster = unique_clusters,
    Num_Expressed_Genes = genes_per_cluster,
    row.names = NULL
)
message("Number of expressed genes after pseudobulk:  ")
gene_counts_df
#               Cluster Num_Expressed_Genes
# 1 Sp03D01                8803
# 2 Sp03D02                8803
# 3 Sp03D03                8803


message('Pseudobulk completed ', k_nice)
message(
  "Dimensions of summed data: ",
  paste(dim(spe_pseudo_k), collapse = " x ")
)

# Rename ncells to nspots

colData(spe_pseudo_k)$nspots <- colData(spe_pseudo_k)$ncells
colData(spe_pseudo_k)$ncells <- NULL # Remove the old column
# add sample-ids to pseudobulk object
colnames(spe_pseudo_k) <- spe_pseudo_k$sample_id


# # Exploring nspots
# 
# min_nspots <- 10
# message("Total nspots: ", sum(spe_pseudo_k$nspots))
# message(
#   "Number of groups with nspots < ",
#   min_nspots,
#   ": ",
#   sum(spe_pseudo_k$nspots < min_nspots)
# )
# message("Summary of nspots:")
# print(summary(spe_pseudo_k$nspots))
# 
# ## Adapted from https://github.com/LieberInstitute/spatialLIBD/blob/devel/R/registration_pseudobulk.R#L137-L154
# ## Drop pseudo-bulked samples that had low initial contribution of raw-samples.
# ## That is, pseudo-bulked samples that are not benefiting from the pseudo-bulking process to obtain higher counts.
# if (!is.null(min_nspots)) {
#   message(
#     Sys.time(),
#     " dropping ",
#     sum(spe_pseudo_k$nspots < min_nspots),
#     " pseudo-bulked samples that are below 'min_nspots'."
#   )
#   spe_pseudo_k <- spe_pseudo_k[, spe_pseudo_k$nspots >= min_nspots]
# }

# Compute mitochondrial expression ratio & other meta-data

is_mito <- which(seqnames(spe_pseudo_k) == "chrM")
spe_pseudo_k$expr_chrM <- colSums(counts(spe_pseudo_k)[is_mito, , drop = FALSE])
spe_pseudo_k$sum_umi <- colSums(counts(spe_pseudo_k))
spe_pseudo_k$expr_chrM_ratio <- spe_pseudo_k$expr_chrM / spe_pseudo_k$sum_umi
spe_pseudo_k$age <- as.numeric(spe_pseudo_k$age)

# Convert other relevant variables to factors

spe_pseudo_k$brain_area2 <- factor(spe_pseudo_k$brain_area2)
if (is.factor(spe_pseudo_k$BayesSpace_pseudo)) {
  spe_pseudo_k$brain_area2 <- droplevels(spe_pseudo_k$brain_area2)
}
levels(spe_pseudo_k$brain_area2)

## Compute the logcounts

message(Sys.time(), " Normalize data and compute PCA")

assays(spe_pseudo_k)

logcounts(spe_pseudo_k) <-
  edgeR::cpm(edgeR::calcNormFactors(spe_pseudo_k), log = TRUE, prior.count = 1)

dim(reducedDim(spe_pseudo_k))

## prepare table of cell counts per sample_id and cluster

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

message(' Normalization completed ')


## Simplify the colData()  for the pseudo-bulked data

# colnames(colData(spe_pseudo_k))
#spe_pseudo_k$BayesSpace
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

## Compute some reduced dims
message('Processing MDS and scarter runPCA')

# keep ncomponents below the number of cells
ncomponents = min(50, ncol(spe_pseudo_k) - 1, nrow(spe_pseudo_k) - 1)

spe_pseudo_k <- scater::runMDS(
  spe_pseudo_k,
  name = "runMDS",
  ncomponents = ncomponents
)


# Fix error when computing runPCA with irlba  - scater default -  too many components relative to the number of features
# -- force exact SVD passing a BiocSingularParam object

spe_pseudo_k <- scater::runPCA(
  spe_pseudo_k,
  ncomponents = ncomponents,
  BSPARAM = ExactParam()
)

dim(reducedDim(spe_pseudo_k, "PCA"))
# [1] 24 23


## For the spatialLIBD shiny app

rowData(spe_pseudo_k)$gene_search <-
  paste0(
    rowData(spe_pseudo_k)$gene_name,
    "; ",
    rowData(spe_pseudo_k)$gene_id
  )

message(
  ' Saving pseudobulk with reduced dims for BS k=',
  k,
  " using brain_area as var for registration"
)

## save RDS file
saveRDS(
  spe_pseudo_k,
  file = file.path(
    dir_rdata,
    paste0("sce_pseudo_PCA_brain_area_k", k_nice, ".rds")
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

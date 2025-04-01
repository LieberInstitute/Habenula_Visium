#################### Prepare enrichment t-stats from Multiome-snRNAseq ##########################################

library("here")
library("Seurat")
library("SingleCellExperiment")
library("spatialLIBD")
library("sessioninfo")

## set hard path to Habenula multiome project WNN Ledien knn=30 resolution=2
dir_outRDS <- here("processed-data", "05_snRNA-seq_model_stats")
inputRDS <- "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/05_Clustering_ARCr/05_rename_idents"
outputRDS <- "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/05_Clustering_ARCr/08_wnn_gene_expression_plts_renamed_idents/"

if (!dir.exists(outputRDS)) {
  dir.create(outputRDS, showWarnings = FALSE, recursive = TRUE)
}

## Read seurat object
rds_name <- here(
  inputRDS,
  "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2.rds"
)
SeuratOBJ <- readRDS(rds_name)
# An object of class Seurat
# 299552 features across 55702 samples within 2 assays
# Active assay: RNA (36601 features, 2000 variable features)
# 3 layers present: data, counts, scale.data
# 1 other assay present: ATAC
# 13 dimensional reductions calculated: pca, umap.unintegrated, integrated.cca, umap, integrated.harmony, lsi, umap.lsi.unintegrated, umap.integrated, tsne.integrated, integrated.lsi.harmony, umap.lsi.integrated, tsne.lsi.integrated, wnn.umap

message("WNN clustering loaded!\nCells: ", length(Cells(x = SeuratOBJ)))
message(
  "WNN containing ",
  nrow(unique(SeuratOBJ[["seurat_clusters"]])),
  " clusters"
)


## Identify and subset clusters annotated as putative `habenula`. Use length of cluster ID as criteria

# ## extract clusters IDs
#
# message("Cluster-IDs from `WNN`")
#
# SeuOBJ_clusters <- Idents(SeuratOBJ)
# head(SeuOBJ_clusters)
# hb_clusters <- unlist(levels(SeuOBJ_clusters))
# hb_clusters
#
# ## filter hb clusters only
#
# no_hb_clust = list()
# for (idx in seq_along(hb_clusters)) { if (nchar(hb_clusters[idx]) <= 4) { no_hb_clust <- append(no_hb_clust, hb_clusters[idx]) } }
# no_hb_clust <- c(unlist(no_hb_clust))
#
# ## Additionally, I make a manual selection of Hb clusters with low-Hb to be removed
# #   - based on % of Hb cells contained in the clusters. More details: https://github.com/LieberInstitute/Hb_multiome/blob/0275ce2f6824b8f22a6efcb1acc9543ca4e1f195/data/06_FULL_SUMMARY_LEIDENr2_knn30_All-Cluster_All-CellTypes_02_11_2025_v3.pdf
#
## join all the cluster with no-hab or low-hab signal
# tmp_hb <- hb_clusters[! hb_clusters %in% c(no_hb_clust)]
# low_hb_clusters <- tmp_hb[c(1, 2, 5, 8, 11, 16)]
# no_hb_clust <- append(no_hb_clust, low_hb_clusters)
#
# hb_clusters <- hb_clusters[! hb_clusters %in% c(no_hb_clust)]
# as.vector(hb_clusters)
# # [1] "C.05 DD_LHb (4.97%)" "C.07 DD_MHb (4.69%)" "C.10 DD_MHb (4.21%)"
# # [4] "C.11 DD_MHb (3.97%)" "C.14 DD_MHb (3.64%)" "C.16 DD_MHb (2.88%)"
# # [7] "C.18 DD_LHb (2.76%)" "C.23 DD_LHb (2.28%)" "C.24 DD_LHb (1.48%)"
# # [10] ...
# length(hb_clusters)
#
# ## subset clusters with Hb cell-types
#
# unique(Idents(SeuratOBJ))
# SeuratOBJ <- subset(SeuratOBJ, idents = as.vector(hb_clusters))
# levels(SeuratOBJ)
# length(Cells(SeuratOBJ))
# # [1] 17827
#
#
# ## removed unused levels
#
# colnames(SeuratOBJ@meta.data)
# head(SeuratOBJ$seurat_clusters)
# SeuratOBJ@meta.data$seurat_clusters <- droplevels(SeuratOBJ@meta.data$seurat_clusters)
#
# # check successful
# levels(SeuratOBJ)
# levels(SeuratOBJ@meta.data$seurat_clusters)
# # [1] "5"  "7"  "10" "11" "14" "16" "18" "23" "24" "30" "33" "36" "40"
#
#

## rename levels in a fancy way to easy Hb clusterID identification

## Faster mode to rename all cluster idents (shorter the name)

levels(SeuratOBJ)
# [1] "C.01 DD_LHb (7.01%)"   "C.02"                  "C.03"
# [4] "C.04 DD_LHb (5.23%)"   "C.05 DD_LHb (4.97%)"   "C.06"
# [7] "C.07 DD_MHb (4.69%)"   "C.08 DD_LHb (4.55%)"   "C.09"
# [10] "C.10 DD_MHb (4.21%)"   "C.11 DD_MHb (3.97%)"   "C.12 DD_LHb (3.93%)"
# [13] "C.13"                  "C.14 DD_MHb (3.64%)"   "C.15"
# ...
new_clust_name <- substr(levels(SeuratOBJ), 1, 11)
names(new_clust_name) <- levels(SeuratOBJ)
SeuratOBJ <- RenameIdents(SeuratOBJ, new_clust_name)
newname_clusters <- levels(SeuratOBJ)
# [1] "C.01 DD_LHb" "C.02"        "C.03"        "C.04 DD_LHb" "C.05 DD_LHb"
# [6] "C.06"        "C.07 DD_MHb" "C.08 DD_LHb" "C.09"        "C.10 DD_MHb"
# [11] "C.11 DD_MHb" "C.12 DD_LHb" "C.13"        "C.14 DD_MHb" "C.15"
# [16] "C.16 DD_MHb" "C.17 LB_Hb " "C.18 DD_LHb" "C.19"        "C.20"
# ...

## Hb clusters mixed with other cell-types will be removed

low_hb_clusters <- c(1, 4, 8, 12, 17, 32)
low_hb_clusters <- paste0("C.", str_pad(low_hb_clusters, width = 2, pad = "0"))
# [1] "C.01" "C.04" "C.08" "C.12" "C.17" "C.32"
low_c <- paste0("^", low_hb_clusters[1], "*")
# low_c = "^C\\.01\\s\\w*"   # "C.01 DD_LHb"
# low_c %in% c(newname_clusters)

##  Rename specific identity classes

SeuratOBJ <- RenameIdents(SeuratOBJ, "C.01 DD_LHb" = "C.01")
SeuratOBJ <- RenameIdents(SeuratOBJ, "C.04 DD_LHb" = "C.04")
SeuratOBJ <- RenameIdents(SeuratOBJ, "C.08 DD_LHb" = "C.08")
SeuratOBJ <- RenameIdents(SeuratOBJ, "C.12 DD_LHb" = "C.12")
SeuratOBJ <- RenameIdents(SeuratOBJ, "C.17 LB_Hb " = "C.17")
SeuratOBJ <- RenameIdents(SeuratOBJ, "C.32 DD_LHb" = "C.32")
levels(SeuratOBJ)

## rename all Seurat clusters as idents

SeuratOBJ$seurat_clusters <- Idents(SeuratOBJ)
unique(SeuratOBJ$seurat_clusters)

## re-order clusters

newname_clusters <- levels(SeuratOBJ)

## extract hb and not hb clusters and sort them

no_hb_clust = list()
desired_order = list()
hb_clusters = list()
for (idx in newname_clusters) {
  if (nchar(idx) <= 4) {
    no_hb_clust <- append(no_hb_clust, idx)
  }
}
no_hb_clust <- sort(c(unlist(no_hb_clust)))
hb_clusters <- sort(newname_clusters[!newname_clusters %in% c(no_hb_clust)])
desired_order <- c(hb_clusters, no_hb_clust)

# Set new identity order

Idents(SeuratOBJ) <- factor(Idents(SeuratOBJ), levels = desired_order)
levels(SeuratOBJ)
# [1] "C.05 DD_LHb" "C.07 DD_MHb" "C.10 DD_MHb" "C.11 DD_MHb" "C.14 DD_MHb"
# [6] "C.16 DD_MHb" "C.18 DD_LHb" "C.23 DD_LHb" "C.24 DD_LHb" "C.30 DD_LHb"
# [11] "C.33 DD_LHb" "C.36 DD_MHb" "C.40 DD_LHb" "C.01"        "C.02"
# [16] "C.03"        "C.04"        "C.06"        "C.08"        "C.09"
# [21] "C.12"        "C.13"        "C.15"        "C.17"        "C.19"
# [26] "C.20"        "C.21"        "C.22"        "C.25"        "C.26"
# [31] "C.27"        "C.28"        "C.29"        "C.31"        "C.32"
# [36] "C.34"        "C.35"        "C.37"        "C.38"        "C.39"
# [41] "C.41"        "C.42"

## rename all Seurat clusters as idents

SeuratOBJ$seurat_clusters <- Idents(SeuratOBJ)
unique(SeuratOBJ$seurat_clusters)
## Check chromatin assay wasn't alterated
class(SeuratOBJ[["ATAC"]])

## save RDS
rds_file_name <- paste0(
  outputRDS,
  "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium.rds"
)
saveRDS(SeuratOBJ, rds_file_name)

message("New seurat with clusters renamed saved on Hb multiome project!")


###################### Retrieve Ensembl IDs for Gene Symbols. ######################
## Note.Cell Ranger ARC reanalyze are barcodes identified as valid barcodes from both
#       `cell-ranger-count` (rna) and `cell-ranger-atac` pipelines run separately

## We use the same reference used on cellranger pipelines

# Cellranger-ARC
reference_gtf <- "/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz"

## Read in the gene information from the annotation GTF file

gtf <- rtracklayer::import(reference_gtf)
#class(gtf)
#"GenomicRanges"

gtf <- gtf[gtf$type == "gene"]
head(gtf)
names(gtf) <- gtf$gene_id # ensembl ids
# Here:
# gene_name = gene symbols (multiome)
# gene_id = gene ensembl

# Extract gene symbols from Seurat object

gene_symbols <- rownames(SeuratOBJ) # Modify if needed for different slot
head(gene_symbols)
# [1] "MIR1302-2HG" "FAM138A"     "OR4F5"       "AL627309.1"  "AL627309.3"

length(gene_symbols)
# [1] 36601

## Import RNA assay into sce object

SeuratOBJ <- DietSeurat(SeuratOBJ, dimreducs = NULL)
sce <- as.SingleCellExperiment(SeuratOBJ, assay = "RNA")
rowData(sce)

### match gene symbols to Ensembl IDs in the sce object

rowData(sce)$gene_id <- unname(gtf$gene_id[match(rownames(sce), gtf$gene_name)])
rowData(sce)$gene_symbol <- rownames(sce)
table(is.na(rownames(sce)))
length(rownames(sce)) # [1] 36601
head(rownames(sce))
# [1] "MIR1302-2HG" "FAM138A"     "OR4F5"       "AL627309.1"  "AL627309.3"

## some validations

table(is.na(rowData(sce)$gene_id)) # ensembl
table(is.na(rowData(sce)$gene_symbol))
table(rownames(sce) %in% gtf$gene_name)
# FALSE  TRUE
# 10    36591
## check genes duplicated or genes that does not match the reference
dup_genes <- setdiff(rownames(sce), gtf$gene_name)
# [1] "TBCE.1"           "LINC01238.1"      "CYB561D2.1"       "MATR3.1"
# [5] "LINC01505.1"      "HSPA14.1"         "GOLGA8M.1"        "GGT1.1"
# [9] "ARMCX5-GPRASP2.1" "TMSB15B.1"

## manage duplicated genes

if (length(dup_genes) > 1) {
  ## remove .1 from gene name
  for (gen in dup_genes) {
    gen_new_name = sub('\\.1', '', gen)
    # print(gen_new_name)
    rownames(sce)[rownames(sce) == gen] <- gen_new_name
  }
  ## match reference again
  rowData(sce)$gene_id <- unname(gtf$gene_id[match(
    rownames(sce),
    gtf$gene_name
  )])
  rowData(sce)$gene_symbol <- rownames(sce)
  ## verify
  dup_genes <- setdiff(rownames(sce), gtf$gene_name)
}

if (length(dup_genes) == 0) {
  message("Annotation ready!")
} else {
  stop()
}

## Perform the enrichment t-stats

sce_modeling_results <- registration_wrapper(
  sce = sce,
  var_registration = "seurat_clusters",
  var_sample_id = "orig.ident",
  gene_ensembl = "gene_id", # gene ensembl ids
  gene_name = "gene_symbol" # gene_names
)


## check out table on enrichment t-statistics

colnames(sce_modeling_results$enrichment)
sce_modeling_results$enrichment[1:5, 1:5]
#                   t_stat_C.01.DD_LHb t_stat_C.02 t_stat_C.03 t_stat_C.04.DD_LHb
# ENSG00000238009          1.4486161  -1.5812290   1.9542309          1.9335732
# ENSG00000241860          1.6045860  -1.2039448   1.2864945          2.1393666
# ENSG00000237491          2.6221958  -0.7505787  -0.1465673          1.5780745

## V2 is a 'polished' and 'sorted' reference version. Here some clusters with low hb cells were removed
##        and the dataset was arranged by first the hb-groups followed by the other cell-types

saveRDS(
  sce_modeling_results,
  here(dir_outRDS, "enrichment_snRNA-multiome_v2b.rds")
)

message(" rna-multiome reference completed!")

# library("slurmjobs")
#
# ## A regular job with 10 cores on the 'imaginary' partition
# job_single("02_multiome_rna_reference", cores = 2, partition = "katun", create_shell = TRUE)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

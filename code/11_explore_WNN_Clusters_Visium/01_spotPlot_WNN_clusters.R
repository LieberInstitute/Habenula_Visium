########################################################################
## Plot spotPlots using multi_gene function on the top5 multiome RNA WNN clusters
##
## Authors. CSC
## Date. April 02, 2025
## Recommended resources on interactive mode: srun --pty --mem=80GB --x11 bash
########################################################################

library("spatialLIBD")
library("Seurat")
library("ggplot2")
library("viridisLite")
library("purrr")
library("tidyverse")
library("here")

## input directories

here()

inputRDS_Dir <- here(
  "processed-data",
  "11_explore_WNN_Clusters_Visium"
)
plotDir <- here(
  "plots",
  "11_explore_WNN_Clusters_Visium",
  "WNN_marker_genes_exploratory_top5"
)
inputCVS_Dir <- here(
  "processed-data",
  "11_explore_WNN_Clusters_Visium"
)

## set path to read sce visium object to plot the top deg 

inputSCE_Dir <- here("processed-data", "04_harmony_BayesSpace", "spe_qcED_spatialLIBD_log.rds")

spe <- readRDS(inputSCE_Dir)
class(spe)
unique(spe$sample_id)
## Quick exploration
cat(" Number of spots:", dim(spe)[2], "\n")

## Check directories

if (!dir.exists(plotDir)) {
  dir.create(plotDir)
}

# For inputRDS_Dir, clusters renamed for Spatial-Registration on Visium project
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium.rds"
seurat_name <- here(inputRDS_Dir, Seurat_base_name)

# Load Seurat to extract `Hb` clustering IDs

SeuratOBJ <- readRDS(here(inputRDS_Dir, Seurat_base_name))
DefaultAssay(SeuratOBJ) <- "RNA"
levels(SeuratOBJ)
## Levels should be
# [1] "C.05 DD_LHb" "C.07 DD_MHb" "C.10 DD_MHb" "C.11 DD_MHb" "C.14 DD_MHb"
# [6] "C.16 DD_MHb" "C.18 DD_LHb" "C.23 DD_LHb" "C.24 DD_LHb" "C.30 DD_LHb"
# [11] "C.33 DD_LHb" "C.36 DD_MHb" "C.40 DD_LHb" "C.01"        "C.02"
# [16] "C.03"        "C.04"        "C.06"        "C.08"        "C.09"
# [21] "C.12"        "C.13"        "C.15"        "C.17"        "C.19"
# [26] "C.20"        "C.21"        "C.22"        "C.25"        "C.26"
# [31] "C.27"        "C.28"        "C.29"        "C.31"        "C.32"
# [36] "C.34"        "C.35"        "C.37"        "C.38"        "C.39"
# [41] "C.41"        "C.42"


## Read DEG from WNN to prepare data to make spotPlots of the top 5 genes highly expressed

# All DEG
DEG_file_name <- "WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_cellTypes_integrated_top50.csv"
DEG_file_name <- here(inputCVS_Dir, DEG_file_name)
df_cluster_names <- read.csv(DEG_file_name)
df_cluster_names <- df_cluster_names |> drop_na(cell_type)
head(df_cluster_names)
# p_val avg_log2FC pct.1 pct.2 p_val_adj cluster     gene            cell_type
# 1     0   4.336764 0.938 0.100         0       1 OTX2-AS1        DD_Inhib.Thal
# 2     0   3.849508 0.900 0.083         0       1      KIT        DD_Inhib.Thal
# 3     0   3.682782 0.927 0.122         0       1    MEIS2      LB_Thalamus/MDm

## Subset clusters annotated as putative `habenula`. Use length of cluster ID as criteria
## extract clusters IDs

message("Filtering habenula cluster-IDs")

SeuOBJ_clusters <- Idents(SeuratOBJ)
remove("SeuratOBJ")

hb_clusters <- unlist(levels(SeuOBJ_clusters))

## Get top 5. Filter habenula clusters only

no_hb_clust = list()
for (idx in seq_along(hb_clusters)) {
  if (nchar(hb_clusters[idx]) <= 4) {
    no_hb_clust <- append(no_hb_clust, hb_clusters[idx])
  }
}
no_hb_clust <- c(unlist(no_hb_clust))
hb_clusters <- hb_clusters[!hb_clusters %in% c(no_hb_clust)]
as.vector(hb_clusters)
# [1] "C.05 DD_LHb" "C.07 DD_MHb" "C.10 DD_MHb" "C.11 DD_MHb" "C.14 DD_MHb"
# [6] "C.16 DD_MHb" "C.18 DD_LHb" "C.23 DD_LHb" "C.24 DD_LHb" "C.30 DD_LHb"
# [11] "C.33 DD_LHb" "C.36 DD_MHb" "C.40 DD_LHb"
# length(hb_clusters)
hb_clusters <- as.integer(trimws(substr(hb_clusters, 3, 4)))
# [1]  5  7 10 11 14 16 18 23 24 30 33 36 40
# get top 5
unique(df_cluster_names$cluster)
top5 <- df_cluster_names |>
  filter(cluster %in% hb_clusters) |>
  group_by(cluster) |>
  top_n(n = 5, wt = avg_log2FC) |>
  select(cluster, gene)

dim(top5)
head(top5)
# cluster gene   
#   <int> <chr>  
# 1       5 RFTN1  
# 2       5 CBLN2  
# 3       5 GALR1  
# 4       5 HTR4   
# 5       5 COL25A1
# 6       7 GPR149 


# unique(top5$cluster)

## Initials for manage plots layout

var_height <- 24 # 24/3=8
var_width <- 36 # 36/4=9
var_point_size <- 3.5

print("Plotting spotPlots with multi-gene function for all the `hb` clusters using he top 5 DEG")

for (clus in as.vector(hb_clusters)) {
  # testing: clus = 5
  message("Processing hb cluster: ", clus)

  # create a vector with gene-ids  
  genes_lst <- top5 |>
    filter(cluster == clus)
  genes_lst <- genes_lst$gene
  # Extract Ensembl ID
  lst_genes <- rowData(spe)$gene_search[
    rowData(spe)$gene_name %in% genes_lst
  ]
  ## Our list of genes
  #lst_genes
  # [1] "RFTN1; ENSG00000131378"   "COL25A1; ENSG00000188517"
  # [3] "HTR4; ENSG00000164270"    "CBLN2; ENSG00000141668"  
  # [5] "GALR1; ENSG00000166573"

  ## define multi-gene method to plot and file names for each cluster
  lst_multi_g <- c(
    z_score = paste0("ZScores_", paste0("C", str_pad(clus, width = 2, pad = "0")), "_multi_genes.pdf"),
    pca = paste0("PCA_", paste0("C", str_pad(clus, width = 2, pad = "0")), "_multi_genes.pdf"),
    sparsity = paste0("sparcity_", paste0("C", str_pad(clus, width = 2, pad = "0")), "_multi_genes.pdf")
  )
  
  map2(
    as.vector(names(lst_multi_g)),
    as.vector(lst_multi_g),
   ~ vis_grid_gene(
      spe = spe,
      geneid = lst_genes,
      multi_gene_method = .x, # z-score, pca, sparcity
      height = var_height,
      width = var_width,
      point_size = var_point_size,
      cont_colors = viridisLite::viridis(21, direction = 1),
      pdf = here(plotDir, .y),
      assayname = "logcounts")
  )
  
}



## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
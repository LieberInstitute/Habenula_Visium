library("here")
library("Seurat")
library("SingleCellExperiment")
library("biomaRt")
# library(EnsDb.Hsapiens.v86)
# library(BSgenome.Hsapiens.UCSC.hg38)
library("spatialLIBD")
library("sessioninfo")

## set hard path to Habenula multiome project WNN Ledien knn=30 resolution=2 
inputRDS <- "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/05_Clustering_ARCr/05_rename_idents"

## Read seurat object
rds_name <- here(inputRDS, "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2.rds")
SeuratOBJ <- readRDS(rds_name)
# An object of class Seurat 
# 299552 features across 55702 samples within 2 assays 
# Active assay: RNA (36601 features, 2000 variable features)
# 3 layers present: data, counts, scale.data
# 1 other assay present: ATAC
# 13 dimensional reductions calculated: pca, umap.unintegrated, integrated.cca, umap, integrated.harmony, lsi, umap.lsi.unintegrated, umap.integrated, tsne.integrated, integrated.lsi.harmony, umap.lsi.integrated, tsne.lsi.integrated, wnn.umap

message("WNN clustering loaded!\nCells: ", length(Cells(x = SeuratOBJ)))
message("WNN containing ", nrow(unique(SeuratOBJ[["seurat_clusters"]])), " clusters")

# brief exploration of the annotations stored in metadata
levels(SeuratOBJ)
# [1] "C.01 DD_LHb (7.01%)"   "C.02"                  "C.03"                 
# [4] "C.04 DD_LHb (5.23%)"   "C.05 DD_LHb (4.97%)"   "C.06"                 
# [7] "C.07 DD_MHb (4.69%)"   "C.08 DD_LHb (4.55%)"   "C.09"                 
# [10] "C.10 DD_MHb (4.21%)"   "C.11 DD_MHb (3.97%)"   "C.12 DD_LHb (3.93%)"  
# ...

## Identified and subset clusters annotated as putative `habenula`. Use length of cluster ID as criteria

## extract clusters IDs

message("Cluster-IDs from `WNN`")
SeuOBJ_clusters <- Idents(SeuratOBJ)
head(SeuOBJ_clusters)
hb_clusters <- unlist(levels(SeuOBJ_clusters))
hb_clusters

## filter hb clusters only 

no_hb_clust = list()

for (idx in seq_along(hb_clusters)) { if (nchar(hb_clusters[idx]) <= 4) { no_hb_clust <- append(no_hb_clust, hb_clusters[idx]) } }
no_hb_clust <- c(unlist(no_hb_clust))

## Additionally, I make a manual selection of Hb clusters with low-Hb to be removed
#   - based on % of Hb cells contained in the clusters. More details: https://github.com/LieberInstitute/Hb_multiome/blob/0275ce2f6824b8f22a6efcb1acc9543ca4e1f195/data/06_FULL_SUMMARY_LEIDENr2_knn30_All-Cluster_All-CellTypes_02_11_2025_v3.pdf 

## join all the cluster with no-hab or low-hab signal
tmp_hb <- hb_clusters[! hb_clusters %in% c(no_hb_clust)]
low_hb_clusters <- tmp_hb[c(1, 2, 5, 8, 11, 16)]
no_hb_clust <- append(no_hb_clust, low_hb_clusters)

hb_clusters <- hb_clusters[! hb_clusters %in% c(no_hb_clust)]
as.vector(hb_clusters)
# [1] "C.05 DD_LHb (4.97%)" "C.07 DD_MHb (4.69%)" "C.10 DD_MHb (4.21%)"
# [4] "C.11 DD_MHb (3.97%)" "C.14 DD_MHb (3.64%)" "C.16 DD_MHb (2.88%)"
# [7] "C.18 DD_LHb (2.76%)" "C.23 DD_LHb (2.28%)" "C.24 DD_LHb (1.48%)"
# [10] ...
length(hb_clusters)

## subset clusters with Hb cell-types

unique(Idents(SeuratOBJ))
SeuratOBJ <- subset(SeuratOBJ, idents = as.vector(hb_clusters))
levels(SeuratOBJ)
length(Cells(SeuratOBJ))
# [1] 17827


## removed unused levels

colnames(SeuratOBJ@meta.data)
head(SeuratOBJ$seurat_clusters)
SeuratOBJ@meta.data$seurat_clusters <- droplevels(SeuratOBJ@meta.data$seurat_clusters)

# check successful
levels(SeuratOBJ)
levels(SeuratOBJ@meta.data$seurat_clusters)
# [1] "5"  "7"  "10" "11" "14" "16" "18" "23" "24" "30" "33" "36" "40"


## rename levels in a fancy way to easy Hb clusterID identification

oldname_clusters <- SeuratOBJ@meta.data$seurat_clusters
newname_clusters <- paste0("HbM.C.", oldname_clusters)
head(newname_clusters)
SeuratOBJ@meta.data$seurat_clusters <- newname_clusters
#levels(SeuratOBJ)
head(SeuratOBJ@meta.data$seurat_clusters)
# [1] "HbM.C.11" "HbM.C.23" "HbM.C.11" "HbM.C.5"  "HbM.C.11" "HbM.C.11"


## Retrieve Ensembl IDs for Gene Symbols
## need to be polish, some gene_id(s) does not match the gene-ensembl id(s)

# Connect to Ensembl database

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")  # For human genes

# Extract gene symbols from Seurat object

gene_symbols <- rownames(SeuratOBJ)  # Modify if needed for different slot
head(gene_symbols)
# [1] "MIR1302-2HG" "FAM138A"     "OR4F5"       "AL627309.1"  "AL627309.3" 

length(gene_symbols)
# [1] 36601

# Convert gene symbols to Ensembl IDs

annotations <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",
  values = gene_symbols,
  mart = mart
)

dim(annotations)
# [1] 26664     2
head(annotations)
#   hgnc_symbol ensembl_gene_id
# 1     A3GALT2 ENSG00000184389
# 2     AADACL3 ENSG00000188984
# 3     AADACL4 ENSG00000204518
# 4        AAK1 ENSG00000115977
# 5       ABCA4 ENSG00000198691
# 6      ABCB10 ENSG00000135776

# Merge Ensembl IDs with Seurat object genes
gene_map <- setNames(annotations$ensembl_gene_id, annotations$hgnc_symbol)
length(gene_map)
# [1] 26664
head(unname(gene_map))
# [1] "ENSG00000184389" "ENSG00000188984" "ENSG00000204518" "ENSG00000115977"
# [5] "ENSG00000198691" "ENSG00000135776"



##### identify clusters 

colnames(SeuratOBJ@meta.data)

## samples included 

table(SeuratOBJ[["orig.ident"]])

## number of cells by cluster

table(SeuratOBJ$seurat_clusters)
# HbM.C.10 HbM.C.11 HbM.C.14 HbM.C.16 HbM.C.18 HbM.C.23 HbM.C.24 HbM.C.30 
# 2344     2212     2026     1607     1535     1269      825      213 
# HbM.C.33 HbM.C.36 HbM.C.40  HbM.C.5  HbM.C.7 
# 186      145       84     2771     2610 



## Import RNA assay in sce object

sce <- as.SingleCellExperiment(SeuratOBJ, assay = "RNA")

rowData(sce)

total_unfiltered_cells <- ncol(sce) # cells in cols
total_unfiltered_cells 
# [1] 55702 / 17827
# unname(gene_map[match(rownames(sce), names(gene_map))])

length(gene_map) # [1] 26664
head(rownames(sce))
# [1] "MIR1302-2HG" "FAM138A"     "OR4F5"       "AL627309.1"  "AL627309.3" 

rowData(sce)$gene_id <- unname(gene_map[match(rownames(sce), names(gene_map))])
rowData(sce)$gene_symbol <- rownames(sce)
rowData(sce)
any(is.na(rownames(sce)))
rowData(sce)
#colnames(colData(sce))
# DataFrame with 36601 rows and 2 columns
# gene_id gene_symbol
# <character> <character>
#   MIR1302-2HG ENSG00000243485 MIR1302-2HG
# FAM138A     ENSG00000237613     FAM138A
# OR4F5       ENSG00000186092       OR4F5
# AL627309.1               NA  AL627309.1
# AL627309.3               NA  AL627309.3


## Perform the spatial registration using ensembl genes ( need to be polish)
# sce_modeling_results <- registration_wrapper(
#   sce = sce,
#   var_registration = "seurat_clusters",
#   var_sample_id = "orig.ident",
#   gene_ensembl = "gene_id", # gene ensembl ids
#   gene_name = "gene_symbol" # gene_names 
# )

## Perform the spatial registration using genes symbols from rna-multiome
sce_modeling_results <- registration_wrapper(
  sce = sce,
  var_registration = "seurat_clusters", # Character vector, C1 / t_stat_C12
  var_sample_id = "orig.ident",
  gene_ensembl = "gene_symbol", # gene ensembl ids
  gene_name = "gene_symbol" # gene_names 
)


## check out table on enrichment t-statistics
sce_modeling_results$enrichment[1:5, 1:5]
# t_stat_HbM.C.10 t_stat_HbM.C.11 t_stat_HbM.C.14 t_stat_HbM.C.16
# AL627309.1      -0.1159839      -0.8339642     0.357303106      0.06155168
# AL627309.5       0.3170415      -0.3913761     1.876669973     -0.57031746
# LINC01409        0.5457288      -1.6630007    -1.192579619     -0.01104049
# LINC01128       -1.0905454       0.9819776     0.255499148     -1.86840108
# LINC00115        0.7564815      -0.2602102     0.005601157     -0.90005001





# ## Create output directories
# dir_rdata <- here("processed-data", "05_snRNA-seq_model_stats")
# dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
# 
# ## Load Habenula snRNA-seq data from the
# ## https://github.com/LieberInstitute/Habenula_Pilot
# ## project.
# 
# ## If we were inside that project
# ## https://github.com/LieberInstitute/Habenula_Pilot/blob/85bfdf36505c829523358a98328a3e2b22bf2042/code/99_paper_figs/12_longer_Heatmap.R#L17-L18
# ##
# 
# habenula_pilot_path <- "~/Dropbox/Code/Habenula_Pilot"
# if (!file.exists(habenula_pilot_path)) {
#   habenula_pilot_path <-
#     "/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula"
# }
# input_sce_path <-
#   file.path(
#     habenula_pilot_path,
#     "processed-data",
#     "04_snRNA-seq",
#     "sce_objects",
#     "sce_final.Rdata"
#   )
# load(input_sce_path, verbose = TRUE)
# 
# table(sce_final$final_Annotations)
# # Astrocyte       Endo Excit.Thal Inhib.Thal      LHb.1      LHb.2
# #       538         38       1800       7612        201        266
# #     LHb.3      LHb.4      LHb.5      LHb.6      LHb.7      MHb.1
# #       134        477         83         39       1014        152
# #     MHb.2      MHb.3  Microglia      Oligo        OPC
# #       540         18        145       2178       1796
# 
# sce_final$final_Annotations_broad <-
#   gsub("\\.[0-9]+", "", sce_final$final_Annotations)
# table(sce_final$final_Annotations_broad, useNA = "ifany")
# # Astrocyte       Endo Excit.Thal Inhib.Thal        LHb        MHb
# #       538         38       1800       7612       2214        710
# # Microglia      Oligo        OPC
# #       145       2178       1796
# 
# registration_vars <-
#   c("final_Annotations", "final_Annotations_broad")
# 
# enrichment_tstats <- lapply(registration_vars, function(current_var) {
#   message(Sys.time(), " processing ", current_var)
#   sce_pseudo <-
#     registration_pseudobulk(
#       sce_final,
#       var_registration = current_var,
#       var_sample_id = "RealSample",
#       pseudobulk_rds_file = file.path(dir_rdata, paste0("pseudobulk_", current_var, ".rds"))
#     )
#   
#   registration_mod <-
#     registration_model(sce_pseudo, covars = "Run")
#   
#   block_cor <-
#     registration_block_cor(sce_pseudo, registration_model = registration_mod)
#   
#   results_enrichment <-
#     registration_stats_enrichment(
#       sce_pseudo,
#       block_cor = block_cor,
#       covars = "Run",
#       gene_ensembl = "ID",
#       gene_name = "Symbol"
#     )
#   
#   file_name = file.path(dir_rdata, paste0("enrichment_", current_var, ".rds"))
#   print(paste0('Enrichment model RDS object', file_name))
#   saveRDS(results_enrichment, file = file_name) # file.path(dir_rdata, paste0("enrichment_", current_var, ".rds"))
#   
#   return(results_enrichment)
# })
# #names(enrichment_tstats) <- registration_vars
# colnames(head(enrichment_tstats$final_Annotations_broad))
# head(enrichment_tstats$final_Annotations_broad)
# #                 t_stat_Astrocyte t_stat_Endo t_stat_Excit.Thal
# # ENSG00000238009        0.4941558   0.9298525        1.06714338
# # ENSG00000241860        0.7344461  -2.5144058        1.34551148
# # ENSG00000237491       -0.3227681  -0.1655309        0.66012110
# # ENSG00000228794       -1.3250693  -2.4485440        1.01537446
# # ENSG00000225880        1.1542415   0.2481395        0.02484108
# # ENSG00000230368        0.1713202   0.2497699        1.75613918
# #                 t_stat_Inhib.Thal t_stat_LHb t_stat_MHb t_stat_Microglia
# # ENSG00000238009         1.2926161 -0.6457847  -1.972686       -0.5768099
# # ENSG00000241860         1.7087987 -1.3389169  -1.756039        0.8147626
# # ENSG00000237491         1.5402734  0.5331491  -2.065161       -2.0884211
# # ENSG00000228794         0.9112540  1.0247560  -1.352135        0.3081748
# # ENSG00000225880         0.3622976 -0.7895245  -1.578354        0.6886194
# # ENSG00000230368         1.0299177 -0.9241120  -1.919964       -1.6361879
# #                 t_stat_Oligo  t_stat_OPC p_value_Astrocyte p_value_Endo
# # ENSG00000238009   -0.6875536  0.15240352         0.6238008   0.35782581
# # ENSG00000241860   -1.4088542  1.38287453         0.4668052   0.01589105
# # ENSG00000237491   -0.2446255  1.14172024         0.7484909   0.86933005
# # ENSG00000228794   -0.1951413  0.57213103         0.1923952   0.01865473
# # ENSG00000225880    0.1278139  0.22600230         0.2550025   0.80525040
# # ENSG00000230368    0.9856077 -0.02113169         0.8648043   0.80399772
# #                 p_value_Excit.Thal p_value_Inhib.Thal p_value_LHb
# # ENSG00000238009         0.29206259         0.20328569   0.5219690
# # ENSG00000241860         0.18575023         0.09494879   0.1878762
# # ENSG00000237491         0.51281562         0.13107542   0.5967704
# # ENSG00000228794         0.31579807         0.36741798   0.3114067
# # ENSG00000225880         0.98030053         0.71896584   0.4342874
# # ENSG00000230368         0.08642671         0.30900658   0.3607631
# #                 p_value_MHb p_value_Microglia p_value_Oligo p_value_OPC
# # ENSG00000238009  0.05521193        0.56717844     0.4955587   0.8796056
# # ENSG00000241860  0.08645330        0.41985147     0.1663270   0.1740825
# # ENSG00000237491  0.04518580        0.04292434     0.8079519   0.2601005
# # ENSG00000228794  0.18364849        0.75948825     0.8462346   0.5703126
# # ENSG00000225880  0.12207045        0.49488810     0.8989136   0.8223065
# # ENSG00000230368  0.06174729        0.10935016     0.3300372   0.9832416
# #                 fdr_Astrocyte   fdr_Endo fdr_Excit.Thal fdr_Inhib.Thal
# # ENSG00000238009     0.8637648 0.48210598      0.8088415      0.5970315
# # ENSG00000241860     0.7806029 0.04995283      0.7711246      0.5016208
# # ENSG00000237491     0.9169061 0.90904042      0.8610136      0.5482505
# # ENSG00000228794     0.5555164 0.05595557      0.8154249      0.6731578
# # ENSG00000225880     0.6245792 0.86195194      0.9938470      0.8529653
# # ENSG00000230368     0.9603161 0.86098997      0.6768724      0.6449839
# #                   fdr_LHb   fdr_MHb fdr_Microglia fdr_Oligo   fdr_OPC
# # ENSG00000238009 0.9518003 0.5097506    0.65387831 0.8410937 0.9553117
# # ENSG00000241860 0.9509080 0.5770869    0.52043485 0.6518411 0.8513704
# # ENSG00000237491 0.9555829 0.4820231    0.09619046 0.9452218 0.8513704
# # ENSG00000228794 0.9509080 0.6954191    0.81683188 0.9539368 0.8533436
# # ENSG00000225880 0.9509080 0.6332252    0.58817081 0.9726442 0.9291777
# # ENSG00000230368 0.9509080 0.5241986    0.18930031 0.7680382 0.9942273
# #                 logFC_Astrocyte logFC_Endo logFC_Excit.Thal
# # ENSG00000238009       0.5169832  1.1725697       0.81335916
# # ENSG00000241860       0.8443182 -3.2925692       1.12174363
# # ENSG00000237491      -0.3198813 -0.1993164       0.47987575
# # ENSG00000228794      -1.3037550 -2.7909634       0.74193769
# # ENSG00000225880       1.2065456  0.3197624       0.01942206
# # ENSG00000230368       0.1779030  0.3147986       1.29600976
# #                 logFC_Inhib.Thal  logFC_LHb logFC_MHb logFC_Microglia
# # ENSG00000238009        0.9786517 -0.4940665 -1.573655      -0.6010286
# # ENSG00000241860        1.4062677 -1.1111834 -1.558357       0.9328809
# # ENSG00000237491        1.0940957  0.3864335 -1.551804      -1.9650995
# # ENSG00000228794        0.6670497  0.7450721 -1.056991       0.3082264
# # ENSG00000225880        0.2826911 -0.6098200 -1.294387       0.7250262
# # ENSG00000230368        0.7774936 -0.6963536 -1.519774      -1.6424652
# #                 logFC_Oligo   logFC_OPC         ensembl       gene
# # ENSG00000238009  -0.6272347  0.11714778 ENSG00000238009 AL627309.1
# # ENSG00000241860  -1.3923262  1.14599996 ENSG00000241860 AL627309.5
# # ENSG00000237491  -0.2121026  0.81767947 ENSG00000237491  LINC01409
# # ENSG00000228794  -0.1713219  0.41956817 ENSG00000228794  LINC01128
# # ENSG00000225880   0.1186566  0.17575396 ENSG00000225880  LINC00115
# # ENSG00000230368   0.8849279 -0.01608617 ENSG00000230368     FAM41C
# 
# ## Format matches the input we need for running
# ## https://research.libd.org/spatialLIBD/reference/layer_stat_cor.html
# modeling_results <- fetch_data(type = "modeling_results")
# colnames(head(modeling_results$enrichment))
# head(modeling_results$enrichment)
# #    t_stat_WM t_stat_Layer1 t_stat_Layer2 t_stat_Layer3 t_stat_Layer4
# # 1 -0.6344143    -1.0321320    0.17815008   -0.72835965    1.56703859
# # 2 -2.4758891     1.2232062   -0.87337451    1.93793650    1.33150141
# # 3 -3.0079360    -0.8564572    2.13358520    0.48741121    0.35212807
# # 4 -1.2916584    -0.9494234   -0.94854397    0.56378302   -0.11206713
# # 5  2.3175897     0.6156900    0.11274780   -0.09907566   -0.03376771
# # 6 -2.2686017    -0.6536163   -0.08615251    1.84786166    0.77710957
# #   t_stat_Layer5 t_stat_Layer6  p_value_WM p_value_Layer1 p_value_Layer2
# # 1    -0.2202707     0.7438713 0.527700348      0.3052551     0.85907467
# # 2     0.4773214    -1.6152865 0.015497447      0.2249981     0.38518446
# # 3     0.6071363     0.3832779 0.003557367      0.3944138     0.03607273
# # 4     1.1536114     1.2634726 0.200358650      0.3453884     0.34583092
# # 5    -0.3434730    -2.4629827 0.023143222      0.5399225     0.91052499
# # 6     0.0355838     0.1945122 0.026108214      0.5153146     0.93156963
# #   p_value_Layer3 p_value_Layer4 p_value_Layer5 p_value_Layer6     fdr_WM
# # 1     0.46861049      0.1212213      0.8262455     0.45922628 0.63711486
# # 2     0.05630996      0.1869669      0.6344905     0.11035355 0.03959651
# # 3     0.62735659      0.7257077      0.5455535     0.70257361 0.01107944
# # 4     0.57454588      0.9110629      0.2522414     0.21024494 0.31550730
# # 5     0.92133658      0.9731501      0.7321823     0.01601949 0.05551142
# # 6     0.06847596      0.4394838      0.9717067     0.84628893 0.06148729
# #   fdr_Layer1 fdr_Layer2 fdr_Layer3 fdr_Layer4 fdr_Layer5 fdr_Layer6
# # 1  0.5644497  0.9418694  0.8284720  0.4139767  0.9596911  0.8115244
# # 2  0.4828399  0.6944277  0.3831776  0.4938944  0.9051814  0.5481106
# # 3  0.6380674  0.2127681  0.9022535  0.8830293  0.8698217  0.9150742
# # 4  0.6001298  0.6683198  0.8845131  0.9674045  0.7017599  0.6698502
# # 5  0.7356770  0.9635319  0.9848034  0.9910492  0.9343088  0.2431891
# # 6  0.7192630  0.9733708  0.4139515  0.7028656  0.9936069  0.9633326
# #           ensembl        gene
# # 1 ENSG00000243485 MIR1302-2HG
# # 2 ENSG00000238009  AL627309.1
# # 3 ENSG00000237491  AL669831.5
# # 4 ENSG00000177757      FAM87B
# # 5 ENSG00000225880   LINC00115
# # 6 ENSG00000230368      FAM41C
# 
# ## Reproducibility information
# print("Reproducibility information:")
# Sys.time()
# proc.time()
# options(width = 120)
# session_info()
#################### Prepare enrichment t-stats from Multiome-snRNAseq ##########################################

library("here")
library("Seurat")
library("SingleCellExperiment")
library("rtracklayer")
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
  "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium.rds"
)
SeuratOBJ <- readRDS(rds_name)
SeuratOBJ
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

levels(SeuratOBJ)
# [1] "C.05.DD_LHb"         "C.07.DD_MHb"         "C.10.DD_MHb"        
# [4] "C.11.DD_MHb"         "C.14.DD_MHb"         "C.16.DD_MHb"        
# [7] "C.18.DD_LHb"         "C.23.DD_LHb"         "C.24.DD_LHb"        
# [10] "C.30.DD_LHb"         "C.33.DD_LHb"         "C.36.DD_MHb"        
# [13] "C.40.DD_LHb"         "C.01.undeterminated" "C.02.DD_Oligo"      
# [16] "C.03.undeterminated" "C.04.undeterminated" "C.06.DD_Exit.Thal"  
# [19] "C.08.undeterminated" "C.09.undeterminated" "C.12.undeterminated"
# [22] "C.13.no-match"       "C.15.DD_Exit.Thal"   "C.17.DD_Exit.Thal"  
# [25] "C.19.DD_Inhib.Thal"  "C.20.DD_Astrocyte"   "C.21.DD_Astrocyte"  
# [28] "C.22.undeterminated" "C.25.undeterminated" "C.26.DD_OPC"        
# [31] "C.27.DD_Microglia"   "C.28.DD_Inhib.Thal"  "C.29.DD_Endo"       
# [34] "C.31.DD_Exit.Thal"   "C.32.undeterminated" "C.34.DD_Oligo"      
# [37] "C.35.undeterminated" "C.37.undeterminated" "C.38.DD_Inhib.Thal" 
# [40] "C.39.DD_Inhib.Thal"  "C.41.DD_Microglia"   "C.42.no-match"      


## Assign new cluster_ann column to Seurat
colnames(SeuratOBJ@meta.data)
table(SeuratOBJ[["seurat_clusters"]])
cluster_ann <- as.vector(Idents(SeuratOBJ))
# assign new identities to the Seurat object
SeuratOBJ$cluster_ann <- cluster_ann
table(SeuratOBJ[["cluster_ann"]])


## Import RNA assay into sce object

## verify that RNA assay has gene symbols
head(rownames(SeuratOBJ[["RNA"]]))
# [1] "MIR1302-2HG" "FAM138A"     "OR4F5"       "AL627309.1"  "AL627309.3" 
# [6] "AL627309.2" 

## make slim Seurat with ony RNA modality, and remove reduction to avoid mismatch issues with sce convertion
SeuratOBJ <- DietSeurat(SeuratOBJ, 
                        assays = "RNA",
                        dimreducs = NULL)
SeuratOBJ
Reductions(SeuratOBJ)
# NULL

## convert Seurat object into sce
sce <- as.SingleCellExperiment(SeuratOBJ)
## verification
sce
# class: SingleCellExperiment 
# dim: 36601 55702 
# metadata(0):
#     assays(3): counts logcounts scaledata
# rownames(36601): MIR1302-2HG FAM138A ... AC007325.4 AC007325.2
# rowData names(0):
#     colnames(55702): S04_AAACAGCCAGAATGAC-1 S04_AAACAGCCAGCAAGGC-1 ...
# S09_TTTGTTGGTCATGCAA-1 S09_TTTGTTGGTTGTTCAC-1
# colData names(30): orig.ident nCount_RNA ... merged_cluster ident
# reducedDimNames(0):
#     mainExpName: RNA
# altExpNames(0):

head(rownames(sce))
# [1] "MIR1302-2HG" "FAM138A"     "OR4F5"       "AL627309.1"  "AL627309.3" 
# [6] "AL627309.2" 


###################### Retrieve Ensembl IDs for Gene Symbols. ######################
## Note.Cell Ranger ARC reanalyze are barcodes identified as valid barcodes from both
#       `cell-ranger-count` (rna) and `cell-ranger-atac` pipelines run separately

## We use the same reference used on cellranger pipelines

# Cellranger-ARC reference
reference_gtf <- "/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz"

## Read the gene information from the annotation GTF file
gtf <- rtracklayer::import(reference_gtf)
class(gtf)
#"GenomicRanges"

# Import and filter for genes
gtf <- import(reference_gtf)
gtf <- gtf[gtf$type == "gene"]
head(gtf)

# Ensure gene_name and gene_id are accessible
gene_map <- data.frame(
    gene_id = gtf$gene_id,
    gene_symbol = gtf$gene_name,
    stringsAsFactors = FALSE
)
head(gene_map)
#         gene_id gene_symbol
# 1 ENSG00000243485 MIR1302-2HG
# 2 ENSG00000237613     FAM138A
# 3 ENSG00000186092       OR4F5
# 4 ENSG00000238009  AL627309.1
# 5 ENSG00000239945  AL627309.3
# 6 ENSG00000239906  AL627309.2

# check duplicated gene symbols (keep first)
table(duplicated(gene_map$gene_symbol))
# FALSE  TRUE 
# 36591    10 
gene_map$gene_symbol[duplicated(gene_map$gene_symbol)]
# [1] "TBCE"           "LINC01238"      "CYB561D2"       "MATR3"         
# [5] "LINC01505"      "HSPA14"         "GOLGA8M"        "GGT1"          
# [9] "ARMCX5-GPRASP2" "TMSB15B"  

# remove duplicated gene symbols (keep first match)
gene_map <- gene_map[!duplicated(gene_map$gene_symbol), ]
table(duplicated(gene_map$gene_symbol))
# FALSE 
# 36591 


## Match based on rownames of sce

matched_ids <- gene_map$gene_id[match(rownames(sce), gene_map$gene_symbol)]
head(matched_ids)
matched_symbols <- rownames(sce)  # These are already your gene symbols
head(matched_symbols)

# Assign to rowData
rowData(sce)$gene_symbol <- matched_symbols
rowData(sce)$gene_id <- matched_ids

# Check matched genes, like those having gene.1 suffix
table(is.na(rowData(sce)$gene_id))
# FALSE  TRUE 
# 36591    10 

# Check unmatched genes, likely those having gene.1 suffix
unmatched <- is.na(rowData(sce)$gene_id)
unmatched_genes <- rowData(sce)$gene_symbol[unmatched]
unmatched_genes
# [1] "TBCE.1"           "LINC01238.1"      "CYB561D2.1"       "MATR3.1"         
# [5] "LINC01505.1"      "HSPA14.1"         "GOLGA8M.1"        "GGT1.1"          
# [9] "ARMCX5-GPRASP2.1" "TMSB15B.1"

# Handle unmatched genes, removing those rows from sce
if (length(unmatched_genes > 0)) {
    # Identify genes ending with ".1"
    genes_to_remove <- grep("\\.1$", rownames(sce), value = TRUE)
    
    sce <- sce[!(rownames(sce) %in% genes_to_remove), ]
    
    # Remap using cleaned but original names
    matched_ids <- gene_map$gene_id[match(rownames(sce), gene_map$gene_symbol)]
    matched_symbols <- rownames(sce)
    
    # Assign to rowData
    rowData(sce)$gene_symbol <- matched_symbols
    rowData(sce)$gene_id <- matched_ids
    
}
# Check matched genes, like those having gene.1 suffix
table(is.na(rowData(sce)$gene_id))
# FALSE
# 29690

message("Genes removed ", length(genes_to_remove))
# Genes removed 6911

# View some annotated entries
head(rowData(sce))
#           gene_symbol         gene_id
#           <character>     <character>
# MIR1302-2HG MIR1302-2HG ENSG00000243485
# FAM138A         FAM138A ENSG00000237613
# OR4F5             OR4F5 ENSG00000186092
# AL627309.1   AL627309.1 ENSG00000238009
# AL627309.3   AL627309.3 ENSG00000239945
# AL627309.2   AL627309.2 ENSG00000239906


# names(gtf) <- gtf$gene_id # ensembl ids
# # Here:
# # gene_name = gene symbols (multiome)
# # gene_id = gene ensembl
# 
# # Extract gene symbols from Seurat object
# 
# gene_symbols <- rownames(SeuratOBJ) # Modify if needed for different slot
# head(gene_symbols)
# # [1] "MIR1302-2HG" "FAM138A"     "OR4F5"       "AL627309.1"  "AL627309.3"
# 
# length(gene_symbols)
# # [1] 36601


## match gene symbols to Ensembl IDs in the sce object
# 
# rowData(sce)$gene_id <- unname(gtf$gene_id[match(rownames(sce), gtf$gene_name)])
# rowData(sce)$gene_symbol <- rownames(sce)
# table(is.na(rownames(sce)))
# length(rownames(sce)) # [1] 36601
# head(rownames(sce))
# # [1] "MIR1302-2HG" "FAM138A"     "OR4F5"       "AL627309.1"  "AL627309.3"
# 
# ## some validations
# 
# table(is.na(rowData(sce)$gene_id)) # ensembl
# table(is.na(rowData(sce)$gene_symbol))
# table(rownames(sce) %in% gtf$gene_name)
# # FALSE  TRUE
# # 10    36591
# ## check genes duplicated or genes that does not match the reference
# dup_genes <- setdiff(rownames(sce), gtf$gene_name)
# # [1] "TBCE.1"           "LINC01238.1"      "CYB561D2.1"       "MATR3.1"
# # [5] "LINC01505.1"      "HSPA14.1"         "GOLGA8M.1"        "GGT1.1"
# # [9] "ARMCX5-GPRASP2.1" "TMSB15B.1"
# 
# ## manage duplicated genes
# 
# if (length(dup_genes) > 1) {
#   ## remove .1 from gene name
#   for (gen in dup_genes) {
#     gen_new_name = sub('\\.1', '', gen)
#     # print(gen_new_name)
#     rownames(sce)[rownames(sce) == gen] <- gen_new_name
#   }
#   ## match reference again
#   rowData(sce)$gene_id <- unname(gtf$gene_id[match(
#     rownames(sce),
#     gtf$gene_name
#   )])
#   rowData(sce)$gene_symbol <- rownames(sce)
#   ## verify
#   dup_genes <- setdiff(rownames(sce), gtf$gene_name)
# }
# 
# if (length(dup_genes) == 0) {
#   message("Annotation ready!")
# } else {
#   stop()
# }

## Perform the enrichment t-stats

sce_modeling_results <- registration_wrapper(
  sce = sce,
  #var_registration = "seurat_clusters",
  var_registration = "cluster_ann",
  var_sample_id = "orig.ident",
  gene_ensembl = "gene_id", # gene ensembl ids
  gene_name = "gene_symbol" # gene_names
)
# var_registration "cluster_ann" contains non-syntatic variables: C.13.no-match, C.42.no-match
# converting to C.13.no.match, C.42.no.match 

## check out table on enrichment t-statistics

colnames(sce_modeling_results$enrichment)
sce_modeling_results$enrichment[1:3, 1:5]
#                   t_stat_C.01.undeterminated t_stat_C.02.DD_Oligo
# ENSG00000241860                  1.5666897           -1.2703843
# ENSG00000237491                  2.5705831           -0.8464480
# ENSG00000228794                  1.1042666           -3.4688938
# ENSG00000225880                  0.6325688           -0.3324544
# ENSG00000230368                  0.1943596            0.7143625
#                   t_stat_C.03.undeterminated t_stat_C.04.undeterminated
# ENSG00000241860                  1.2485112                  2.2029260
# ENSG00000237491                 -0.1863302                  1.6434360
# ENSG00000228794                  0.3765422                  1.1600994
# ENSG00000225880                  0.3383681                  0.6374913
# ENSG00000230368                  1.1264133                  2.3185947
#                   t_stat_C.05.DD_LHb
# ENSG00000241860          0.2317114
# ENSG00000237491          0.4301931
# ENSG00000228794          1.1586025
# ENSG00000225880          1.0757814
# ENSG00000230368          0.5298103

## V2 has hb clusters annotated (curated - 13 clusters)
## V3 has all clusters pre-annotated for EDA with Clustering-Registration

saveRDS(
  sce,
  # here(dir_outRDS, "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_v2b.rds")
  here(dir_outRDS, "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_v3.rds")
)

saveRDS(
  sce_modeling_results,
  # here(dir_outRDS, "enrichment_snRNA-multiome_v2b.rds")
  here(dir_outRDS, "enrichment_snRNA-multiome_v3.rds")
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

library("spatialLIBD")
library("tidyverse")
library("here")
library("lobstr")
library("sessioninfo")

source(here("code/04_harmony_BayesSpace", "func_DGE_marker_gene_lists.R")) # Call functions to read paths

### Implementing function to filter the data-driven marker gene list from top50r
# This is de function to call ... get_Top_DataDriven_markers_genes()

s_path_name <- here("data", "sfigu_top_50_MarkerGenes_Table.xlsx")
Hb_gene_markers <- as.data.frame(read_excel(s_path_name, na = "---")) # sheet = "data"
## Note: cellType.target is the column you want to use. That is the "target" cell type that the data corresponds to,
## the second cellType is the second highest non-target cell type (so the cell type we are comparing the target cell type to)

head(Hb_gene_markers)
colnames(Hb_gene_markers)
# [1] "gene"                  "cellType.target"       "mean.target"
# [4] "cellType"              "mean"                  "ratio"
# [7] "rank_ratio"            "Symbol"                "anno_ratio"
# [10] "logFC"                 "log.p.value"           "log.FDR"
# [13] "std.logFC"             "std.logFC_rank_marker" "anno_logFC"

# top50_LHb_genes_byratio <- Hb_gene_markers %>%
#   dplyr::arrange(cellType.target) %>%
#   select(c(cellType.target, rank_ratio, Symbol, log.p.value), 1) %>%   #, starts_with(f)
#   dplyr::filter((cellType.target == 'LHb')) #| (cellType.target == 'MHb') %>%
# #slice_head(n = 50)
#
# top50_MHb_genes_byratio <- Hb_gene_markers %>%
#   dplyr::arrange(cellType.target) %>%
#   select(c(cellType.target, rank_ratio, Symbol, log.p.value), 1) %>%
#   dplyr::filter((cellType.target == 'MHb')) #| (cellType.target == 'MHb') %>%
#
# #nrow(top50_LHb_genes_byratio)
# #tail(top50_LHb_genes_byratio, n=5)
#
# if ( (nrow(top50_LHb_genes_byratio) > 0) &  (nrow(top50_MHb_genes_byratio) > 0) ) {
#   markers.custom = list(
#     'LHb_putative' = top50_LHb_genes_byratio$Symbol, # 'GRIN1','MAP2'),
#     'MHb_putative' = top50_MHb_genes_byratio$Symbol # 'SLC17A8'),
#   )
# }









######################## Load the full marker genes lists (Known broad and Data driven for Habenula) ############################

# We have access to 3 gene markers lists:

# Erik and Top50r putative marker genes merged
markers.custom <- get_erik_and_Hb_markers_genes() # merged lists
# prefix_name <- "all_gm" # prefix to save matched markers found in the clusters
# markers.custom <- get_bukola_markers_genes_Hb()           # Bukola lists
# prefix_name <- 'erik_gm'
# markers.custom <- get_Top50r_markers_genes_Hb()           # Top50r lists (putative Hb)
# prefix_name <- 'Top50r_gm'

# str(markers.custom)
# List of 14
# $ neuron                   : chr [1:2] "SYT1" "SNAP25"
# $ excitatory_neuron        : chr [1:2] "SLC17A6" "SLC17A7"
# $ inhibitory_neuron        : chr [1:2] "GAD1" "GAD2"
# $ mediodorsal thalamus     : chr [1:10] "EPHA4" "PDYN" "LYPD6B" "LYPD6" ...
# $ Hb neuron specific       : chr [1:4] "POU2F2" "POU4F1" "GPR151" "CALB2"
# $ MHB neuron specific      : chr [1:3] "TAC1" "CHAT" "CHRNB4"
# $ LHB neuron specific      : chr [1:2] "HTR2C" "MMRN1"
# $ oligodendrocyte          : chr [1:2] "MOBP" "MBP"
# $ oligodendrocyte_precursor: chr [1:2] "PDGFRA" "VCAN"
# $ microglia                : chr [1:2] "C3" "CSF1R"
# $ astrocyte                : chr [1:2] "GFAP" "AQP4"
# $ Endo/CP                  : chr [1:4] "TTR" "FOLR1" "FLT1" "CLDN5"
# $ MHb_putative             : chr [1:50] "CHAT" "LINC01307" "NEUROD1" "CHRNB4" ...
# $ LHb_putative             : chr [1:50] "HTR4" "BVES" "NRP1" "HTR2C" ...


markers.custom$MHb_putative
# [1] "CHAT"       "LINC01307"  "NEUROD1"    "CHRNB4"     "LINC02143"
# [6] "AC114321.1" "AC104170.1" "AC079760.2" "AC024610.2" "AC022382.2"
# ...

# # markers manually added for testing functions
# new_gm <- c('AQP4', 'MT-ND2')
# markers.custom$MHb <- append(markers.custom$MHb, new_gm)
# markers.custom$MHb

####################### Run multi_gene analysis for exploratory purposes ###########################

## Set directory for data and plots

dir_plots <- here("plots", "", "04_harmony_BayesSpace")
dir_rdata <- here("processed-data", "04_harmony_BayesSpace")

## set path to read RDS object. In this case I set the QCed data with log normalized counts.
spe_in_path <- here("processed-data", "04_harmony_BayesSpace", "spe_qcED_spatialLIBD_log.rds") 
spe <- readRDS(spe_in_path)
unique(spe$sample_id)

## Quick exploration
cat(" Number of spots:", dim(spe)[2], "\n")

## Set some initials for manage plots
var_height <- 24 # 24/3=8
var_width <- 36 # 36/4=9
var_point_size <- 3.5


#######################  Inspect WM known gene markers   #######################

lst_white_matter_genes <- c("GFAP", "AQP4", "MBP", "PLP1")

# Extract Ensembl ID
lst_WM <- rowData(spe)$gene_search[
    rowData(spe)$gene_name %in% lst_white_matter_genes
]

## Our list of white matter genes
lst_WM
# [1] "GFAP; ENSG00000131095" "AQP4; ENSG00000171885" "MBP; ENSG00000197971"
# [4] "PLP1; ENSG00000123560"

## -----------------------------
## in-tissue metrics
# set the number of top genes selected for the multi gene analysis
suffix_name <- ""

## define multi-gene method to plot
lst_multi_g <- c(
    z_score = paste0("literature_multi_genes_Zs_WM", suffix_name, ".pdf"),
    pca = paste0("literature_multi_genes_PCA_WM", suffix_name, ".pdf"),
    sparsity = paste0("literature_multi_genes_Sp_WM", suffix_name, ".pdf")
)

print("Ploting multi-genes for WM gene markers")

map2(as.vector(names(lst_multi_g)), as.vector(lst_multi_g), ~ vis_grid_gene(
    spe = spe,
    geneid = lst_WM,
    multi_gene_method = .x, # z-score, pca, sparcity
    height = var_height,
    width = var_width,
    point_size = var_point_size,
    # cont_colors = viridisLite::turbo(21, direction = 1),
    cont_colors = viridisLite::viridis(21, direction = 1),
    # return_plots = TRUE,
    pdf = here(dir_plots, "Marker_genes_exploratory_initial", .y),
    assayname = "logcounts" #"counts"
))




#######################  Inspect `Hb neuron specific` genes marker #######################

## -----------------------------
# markers.custom$`Hb neuron specific`

suffix_name <- "Hb_neuron.pdf"

# Extract Ensembl ID
lst_Habenula <- rowData(spe)$gene_search[
    rowData(spe)$gene_name %in% markers.custom$`Hb neuron specific`
]

lst_multi_g <- c(
    z_score = paste0("literature_multi_genes_Zs_", suffix_name),
    pca = paste0("literature_multi_genes_PCA_", suffix_name),
    sparsity = paste0("literature_multi_genes_Sp_", suffix_name)
)

print("Ploting multi-gene model for habenula neuron specific gene markers")

map2(as.vector(names(lst_multi_g)), as.vector(lst_multi_g), ~ vis_grid_gene(
    spe = spe,
    geneid = lst_Habenula,
    multi_gene_method = .x,
    height = var_height,
    width = var_width,
    point_size = var_point_size,
    cont_colors = viridisLite::viridis(21, direction = 1),
    pdf = here(dir_plots, "Marker_genes_exploratory_initial", .y),
    assayname = "logcounts"
))



#######################  Inspect `MHB neuron specific` genes marker #######################

## -----------------------------
# markers.custom$`MHB neuron specific` / `LHB neuron specific` / `mediodorsal thalamus`

suffix_name <- "Medial_Hb_neuron.pdf"

# Extract Ensembl ID
lst_Habenula <- rowData(spe)$gene_search[
    rowData(spe)$gene_name %in% markers.custom$`MHB neuron specific`
]

## Multi-gene function crashes if searching have no expression variation (CHAT and CHRNB4) in some of the selected genes
## Thus, we adjusted the searching for only the TAC gene
# lst_Habenula <- lst_Habenula[1:2]
suffix_name <- "Medial_Hb_TAC_neuron.pdf"

lst_multi_g <- c(
    z_score = paste0("literature_multi_genes_Zs_", suffix_name),
    pca = paste0("literature_multi_genes_PCA_", suffix_name),
    sparsity = paste0("literature_multi_genes_Sp_", suffix_name)
)

print("Ploting multi-gene model for Medial Habenula neuron specific gene markers")

map2(as.vector(names(lst_multi_g)), as.vector(lst_multi_g), ~ vis_grid_gene(
    spe = spe,
    geneid = lst_Habenula,
    multi_gene_method = .x,
    height = var_height,
    width = var_width,
    point_size = var_point_size,
    cont_colors = viridisLite::viridis(21, direction = 1),
    pdf = here(dir_plots, "Marker_genes_exploratory_initial", .y),
    assayname = "logcounts"
))



#######################  Inspect `LHB neuron specific` genes marker #######################

## -----------------------------
# markers.custom$`LHB neuron specific` / `mediodorsal thalamus`

suffix_name <- "Lateral_Hb_neuron.pdf"

# Extract Ensembl ID
lst_Habenula <- rowData(spe)$gene_search[
    rowData(spe)$gene_name %in% markers.custom$`LHB neuron specific`
]

lst_multi_g <- c(
    z_score = paste0("literature_multi_genes_Zs_", suffix_name),
    pca = paste0("literature_multi_genes_PCA_", suffix_name),
    sparsity = paste0("literature_multi_genes_Sp_", suffix_name)
)

print("Ploting multi-gene model for Lateral Habenula neuron specific gene markers")

map2(as.vector(names(lst_multi_g)), as.vector(lst_multi_g), ~ vis_grid_gene(
    spe = spe,
    geneid = lst_Habenula,
    multi_gene_method = .x,
    height = var_height,
    width = var_width,
    point_size = var_point_size,
    cont_colors = viridisLite::viridis(21, direction = 1),
    pdf = here(dir_plots, "Marker_genes_exploratory_initial", .y),
    assayname = "logcounts"
))



#######################  Inspect Lateral Habenula genes marker from top50r #######################

## -----------------------------
## ALL 50 genes markers
# markers.custom$LHb_putative

suffix_name <- "LH_ALL.pdf"

# Extract Ensembl ID
lst_LH <- rowData(spe)$gene_search[
    rowData(spe)$gene_name %in% markers.custom$LHb_putative
]

lst_multi_g <- c(
    z_score = paste0("multi_genes_Zs_", suffix_name),
    pca = paste0("multi_genes_PCA_", suffix_name),
    sparsity = paste0("multi_genes_Sp_", suffix_name)
)

print("Ploting multi-gene model for ALL Lateral Habenula gene markers")

map2(as.vector(names(lst_multi_g)), as.vector(lst_multi_g), ~ vis_grid_gene(
    spe = spe,
    geneid = lst_LH,
    multi_gene_method = .x,
    height = var_height,
    width = var_width,
    point_size = var_point_size,
    # cont_colors = viridisLite::turbo(21, direction = 1),
    cont_colors = viridisLite::viridis(21, direction = 1),
    # return_plots = TRUE,
    pdf = here(dir_plots, .y),
    assayname = "counts"
))


## -----------------------------
## Top 25 genes markers

suffix_name <- "LH_Top25.pdf"

# markers.custom$LHb_putative
# Extract Ensembl ID
lst_LH <- rowData(spe)$gene_search[
    rowData(spe)$gene_name %in% markers.custom$LHb_putative[1:25]
]

lst_multi_g <- c(
    z_score = paste0("multi_genes_Zs_", suffix_name),
    pca = paste0("multi_genes_PCA_", suffix_name),
    sparsity = paste0("multi_genes_Sp_", suffix_name)
)

print("Ploting multi-gene model for the top 25 Lateral Habenula gene markers")

map2(as.vector(names(lst_multi_g)), as.vector(lst_multi_g), ~ vis_grid_gene(
    spe = spe,
    geneid = lst_LH,
    multi_gene_method = .x,
    height = var_height,
    width = var_width,
    point_size = var_point_size,
    # cont_colors = viridisLite::turbo(21, direction = 1),
    cont_colors = viridisLite::viridis(21, direction = 1),
    # return_plots = TRUE,
    pdf = here(dir_plots, .y),
    assayname = "counts"
))


#######################  Inspect Medial Habenula gene markers from top 50r #######################

## -----------------------------
## ALL 50 genes markers
# markers.custom$MHb_putative

suffix_name <- "MH_ALL.pdf"

# Extract Ensembl ID
lst_MH <- rowData(spe)$gene_search[
    rowData(spe)$gene_name %in% markers.custom$MHb_putative
]

lst_multi_g <- c(
    z_score = paste0("multi_genes_Zs_", suffix_name),
    pca = paste0("multi_genes_PCA_", suffix_name),
    sparsity = paste0("multi_genes_Sp_", suffix_name)
)

print("Ploting multi-gene model for ALL Medial Habenula gene markers")

map2(as.vector(names(lst_multi_g)), as.vector(lst_multi_g), ~ vis_grid_gene(
    spe = spe,
    geneid = lst_MH,
    multi_gene_method = .x,
    height = var_height,
    width = var_width,
    point_size = var_point_size,
    # cont_colors = viridisLite::turbo(21, direction = 1),
    cont_colors = viridisLite::viridis(21, direction = 1),
    # return_plots = TRUE,
    pdf = here(dir_plots, .y),
    assayname = "counts"
))

## -----------------------------
## Top 25 genes markers

suffix_name <- "MH_Top25.pdf"

# Extract Ensembl ID
lst_MH <- rowData(spe)$gene_search[
    rowData(spe)$gene_name %in% markers.custom$MHb_putative[1:25]
]

lst_multi_g <- c(
    z_score = paste0("multi_genes_Zs_", suffix_name),
    pca = paste0("multi_genes_PCA_", suffix_name),
    sparsity = paste0("multi_genes_Sp_", suffix_name)
)

print("Ploting multi-gene model for the top 25 Medial Habenula gene markers")

map2(as.vector(names(lst_multi_g)), as.vector(lst_multi_g), ~ vis_grid_gene(
    spe = spe,
    geneid = lst_MH,
    multi_gene_method = .x,
    height = var_height,
    width = var_width,
    point_size = var_point_size,
    # cont_colors = viridisLite::turbo(21, direction = 1),
    cont_colors = viridisLite::viridis(21, direction = 1),
    # return_plots = TRUE,
    pdf = here(dir_plots, .y),
    assayname = "counts"
))



#######################  Inspect `mediodorsal thalamus` genes marker #######################

## -----------------------------
# markers.custom$`mediodorsal thalamus`

suffix_name <- "thalamus.pdf"

# Extract Ensembl ID
lst_thalamus <- rowData(spe)$gene_search[
    rowData(spe)$gene_name %in% markers.custom$`mediodorsal thalamus`
]

lst_multi_g <- c(
    z_score = paste0("multi_genes_Zs_", suffix_name),
    pca = paste0("multi_genes_PCA_", suffix_name),
    sparsity = paste0("multi_genes_Sp_", suffix_name)
)

print("Ploting multi-gene model for thalamus gene markers")

map2(as.vector(names(lst_multi_g)), as.vector(lst_multi_g), ~ vis_grid_gene(
    spe = spe,
    geneid = lst_thalamus,
    multi_gene_method = .x,
    height = var_height,
    width = var_width,
    point_size = var_point_size,
    # cont_colors = viridisLite::turbo(21, direction = 1),
    cont_colors = viridisLite::viridis(21, direction = 1),
    # return_plots = TRUE,
    pdf = here(dir_plots, .y),
    assayname = "counts"
))


# ==============================================================================

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


# [1] "Reproducibility information:"
# > Sys.time()
# [1] "2024-07-11 12:09:15 EDT"
# > proc.time()
# user   system  elapsed 
# 68.417    2.262 2990.522 
# > options(width = 120)
# > session_info()
# .8       2022-06-12 [2] CRAN (R 4.3.2)
# benchmarkmeData          1.0.4       2020-04-23 [2] CRAN (R 4.3.2)
# Biobase                * 2.62.0      2023-10-24 [2] Bioconductor
# BiocFileCache            2.10.1      2023-10-26 [2] Bioconductor
# BiocGenerics           * 0.48.1      2023-11-01 [2] Bioconductor
# BiocIO                   1.12.0      2023-10-24 [2] Bioconductor
# BiocManager              1.30.22     2023-08-08 [2] CRAN (R 4.3.2)
# BiocNeighbors            1.20.2      2024-01-07 [2] Bioconductor 3.18 (R 4.3.2)
# BiocParallel             1.36.0      2023-10-24 [2] Bioconductor
# BiocSingular             1.18.0      2023-10-24 [2] Bioconductor
# BiocVersion              3.18.1      2023-11-15 [2] Bioconductor
# Biostrings               2.70.2      2024-01-28 [2] Bioconductor 3.18 (R 4.3.2)
# bit                      4.0.5       2022-11-15 [2] CRAN (R 4.3.2)
# bit64                    4.0.5       2020-08-30 [2] CRAN (R 4.3.2)
# bitops                   1.0-7       2021-04-24 [2] CRAN (R 4.3.2)
# blob                     1.2.4       2023-03-17 [2] CRAN (R 4.3.2)
# bslib                    0.6.1       2023-11-28 [2] CRAN (R 4.3.2)
# cachem                   1.0.8       2023-05-01 [2] CRAN (R 4.3.2)
# cellranger               1.1.0       2016-07-27 [2] CRAN (R 4.3.2)
# cli                      3.6.2       2023-12-11 [2] CRAN (R 4.3.2)
# codetools                0.2-19      2023-02-01 [3] CRAN (R 4.3.2)
# colorspace               2.1-0       2023-01-23 [2] CRAN (R 4.3.2)
# config                   0.3.2       2023-08-30 [2] CRAN (R 4.3.2)
# cowplot                  1.1.3       2024-01-22 [2] CRAN (R 4.3.2)
# crayon                   1.5.2       2022-09-29 [2] CRAN (R 4.3.2)
# curl                     5.2.0       2023-12-08 [2] CRAN (R 4.3.2)
# data.table               1.15.0      2024-01-30 [2] CRAN (R 4.3.2)
# DBI                      1.2.1       2024-01-12 [2] CRAN (R 4.3.2)
# dbplyr                   2.4.0       2023-10-26 [2] CRAN (R 4.3.2)
# DelayedArray           * 0.28.0      2023-10-24 [2] Bioconductor
# DelayedMatrixStats       1.24.0      2023-10-24 [2] Bioconductor
# digest                   0.6.34      2024-01-11 [2] CRAN (R 4.3.2)
# doParallel               1.0.17      2022-02-07 [2] CRAN (R 4.3.2)
# dotCall64                1.1-1       2023-11-28 [2] CRAN (R 4.3.2)
# dplyr                  * 1.1.4       2023-11-17 [2] CRAN (R 4.3.2)
# DT                       0.31        2023-12-09 [2] CRAN (R 4.3.2)
# edgeR                    4.0.14      2024-01-29 [2] Bioconductor 3.18 (R 4.3.2)
# ellipsis                 0.3.2       2021-04-29 [2] CRAN (R 4.3.2)
# ExperimentHub            2.10.0      2023-10-24 [2] Bioconductor
# fansi                    1.0.6       2023-12-08 [2] CRAN (R 4.3.2)
# farver                   2.1.1       2022-07-06 [2] CRAN (R 4.3.2)
# fastmap                  1.1.1       2023-02-24 [2] CRAN (R 4.3.2)
# fields                   15.2        2023-08-17 [2] CRAN (R 4.3.2)
# filelock                 1.0.3       2023-12-11 [2] CRAN (R 4.3.2)
# forcats                * 1.0.0       2023-01-29 [2] CRAN (R 4.3.2)
# foreach                  1.5.2       2022-02-02 [2] CRAN (R 4.3.2)
# generics                 0.1.3       2022-07-05 [2] CRAN (R 4.3.2)
# GenomeInfoDb           * 1.38.5      2023-12-28 [2] Bioconductor 3.18 (R 4.3.2)
# GenomeInfoDbData         1.2.11      2024-02-09 [2] Bioconductor
# GenomicAlignments        1.38.2      2024-01-16 [2] Bioconductor 3.18 (R 4.3.2)
# GenomicRanges          * 1.54.1      2023-10-29 [2] Bioconductor
# ggbeeswarm               0.7.2       2023-04-29 [2] CRAN (R 4.3.2)
# ggplot2                * 3.4.4       2023-10-12 [2] CRAN (R 4.3.2)
# ggrepel                  0.9.5       2024-01-10 [2] CRAN (R 4.3.2)
# glue                     1.7.0       2024-01-09 [2] CRAN (R 4.3.2)
# golem                    0.4.1       2023-06-05 [2] CRAN (R 4.3.2)
# gridExtra                2.3         2017-09-09 [2] CRAN (R 4.3.2)
# gtable                   0.3.4       2023-08-21 [2] CRAN (R 4.3.2)
# HDF5Array              * 1.30.0      2023-10-24 [2] Bioconductor
# here                   * 1.0.1       2020-12-13 [2] CRAN (R 4.3.2)
# hms                      1.1.3       2023-03-21 [2] CRAN (R 4.3.2)
# htmltools                0.5.7       2023-11-03 [2] CRAN (R 4.3.2)
# htmlwidgets              1.6.4       2023-12-06 [2] CRAN (R 4.3.2)
# httpuv                   1.6.14      2024-01-26 [2] CRAN (R 4.3.2)
# httr                     1.4.7       2023-08-15 [2] CRAN (R 4.3.2)
# interactiveDisplayBase   1.40.0      2023-10-24 [2] Bioconductor
# IRanges                * 2.36.0      2023-10-24 [2] Bioconductor
# irlba                    2.3.5.1     2022-10-03 [2] CRAN (R 4.3.2)
# iterators                1.0.14      2022-02-05 [2] CRAN (R 4.3.2)
# jquerylib                0.1.4       2021-04-26 [2] CRAN (R 4.3.2)
# jsonlite                 1.8.8       2023-12-04 [2] CRAN (R 4.3.2)
# KEGGREST                 1.42.0      2023-10-24 [2] Bioconductor
# labeling                 0.4.3       2023-08-29 [2] CRAN (R 4.3.2)
# later                    1.3.2       2023-12-06 [2] CRAN (R 4.3.2)
# lattice                  0.22-5      2023-10-24 [3] CRAN (R 4.3.2)
# lazyeval                 0.2.2       2019-03-15 [2] CRAN (R 4.3.2)
# lifecycle                1.0.4       2023-11-07 [2] CRAN (R 4.3.2)
# limma                    3.58.1      2023-10-31 [2] Bioconductor
# lobstr                 * 1.1.2       2022-06-22 [2] CRAN (R 4.3.2)
# locfit                   1.5-9.8     2023-06-11 [2] CRAN (R 4.3.2)
# lubridate              * 1.9.3       2023-09-27 [2] CRAN (R 4.3.2)
# magick                   2.8.2       2023-12-20 [2] CRAN (R 4.3.2)
# magrittr                 2.0.3       2022-03-30 [2] CRAN (R 4.3.2)
# maps                     3.4.2       2023-12-15 [2] CRAN (R 4.3.2)
# Matrix                 * 1.6-5       2024-01-11 [3] CRAN (R 4.3.2)
# MatrixGenerics         * 1.14.0      2023-10-24 [2] Bioconductor
# matrixStats            * 1.2.0       2023-12-11 [2] CRAN (R 4.3.2)
# memoise                  2.0.1       2021-11-26 [2] CRAN (R 4.3.2)
# mime                     0.12        2021-09-28 [2] CRAN (R 4.3.2)
# munsell                  0.5.0       2018-06-12 [2] CRAN (R 4.3.2)
# paletteer                1.6.0       2024-01-21 [2] CRAN (R 4.3.2)
# pillar                   1.9.0       2023-03-22 [2] CRAN (R 4.3.2)
# pkgconfig                2.0.3       2019-09-22 [2] CRAN (R 4.3.2)
# plotly                   4.10.4      2024-01-13 [2] CRAN (R 4.3.2)
# png                      0.1-8       2022-11-29 [2] CRAN (R 4.3.2)
# promises                 1.2.1       2023-08-10 [2] CRAN (R 4.3.2)
# purrr                  * 1.0.2       2023-08-10 [2] CRAN (R 4.3.2)
# R6                       2.5.1       2021-08-19 [2] CRAN (R 4.3.2)
# rappdirs                 0.3.3       2021-01-31 [2] CRAN (R 4.3.2)
# RColorBrewer             1.1-3       2022-04-03 [2] CRAN (R 4.3.2)
# Rcpp                     1.0.12      2024-01-09 [2] CRAN (R 4.3.2)
# RCurl                    1.98-1.14   2024-01-09 [2] CRAN (R 4.3.2)
# readr                  * 2.1.5       2024-01-10 [2] CRAN (R 4.3.2)
# readxl                 * 1.4.3       2023-07-06 [2] CRAN (R 4.3.2)
# rematch2                 2.1.2       2020-05-01 [2] CRAN (R 4.3.2)
# restfulr                 0.0.15      2022-06-16 [2] CRAN (R 4.3.2)
# rhdf5                  * 2.46.1      2023-11-29 [2] Bioconductor 3.18 (R 4.3.2)
# rhdf5filters             1.14.1      2023-11-06 [2] Bioconductor
# Rhdf5lib                 1.24.2      2024-02-07 [2] Bioconductor 3.18 (R 4.3.2)
# rjson                    0.2.21      2022-01-09 [2] CRAN (R 4.3.2)
# rlang                    1.1.3       2024-01-10 [2] CRAN (R 4.3.2)
# rprojroot                2.0.4       2023-11-05 [2] CRAN (R 4.3.2)
# Rsamtools                2.18.0      2023-10-24 [2] Bioconductor
# RSQLite                  2.3.5       2024-01-21 [2] CRAN (R 4.3.2)
# rsvd                     1.0.5       2021-04-16 [2] CRAN (R 4.3.2)
# rtracklayer              1.62.0      2023-10-24 [2] Bioconductor
# S4Arrays               * 1.2.0       2023-10-24 [2] Bioconductor
# S4Vectors              * 0.40.2      2023-11-23 [2] Bioconductor 3.18 (R 4.3.2)
# sass                     0.4.8       2023-12-06 [2] CRAN (R 4.3.2)
# ScaledMatrix             1.10.0      2023-10-24 [2] Bioconductor
# scales                   1.3.0       2023-11-28 [2] CRAN (R 4.3.2)
# scater                   1.30.1      2023-11-16 [2] Bioconductor
# scuttle                  1.12.0      2023-10-24 [2] Bioconductor
# sessioninfo            * 1.2.2       2021-12-06 [2] CRAN (R 4.3.2)
# shiny                    1.8.0       2023-11-17 [2] CRAN (R 4.3.2)
# shinyWidgets             0.8.1       2024-01-10 [2] CRAN (R 4.3.2)
# SingleCellExperiment   * 1.24.0      2023-10-24 [2] Bioconductor
# spam                     2.10-0      2023-10-23 [2] CRAN (R 4.3.2)
# SparseArray            * 1.2.3       2023-12-25 [2] Bioconductor 3.18 (R 4.3.2)
# sparseMatrixStats        1.14.0      2023-10-24 [2] Bioconductor
# SpatialExperiment      * 1.12.0      2023-10-24 [2] Bioconductor
# spatialLIBD            * 1.15.4      2024-05-01 [1] Github (LieberInstitute/spatialLIBD@77a5303)
# statmod                  1.5.0       2023-01-06 [2] CRAN (R 4.3.2)
# stringi                  1.8.3       2023-12-11 [2] CRAN (R 4.3.2)
# stringr                * 1.5.1       2023-11-14 [2] CRAN (R 4.3.2)
# SummarizedExperiment   * 1.32.0      2023-10-24 [2] Bioconductor
# tibble                 * 3.2.1       2023-03-20 [2] CRAN (R 4.3.2)
# tidyr                  * 1.3.1       2024-01-24 [2] CRAN (R 4.3.2)
# tidyselect               1.2.0       2022-10-10 [2] CRAN (R 4.3.2)
# tidyverse              * 2.0.0       2023-02-22 [2] CRAN (R 4.3.2)
# timechange               0.3.0       2024-01-18 [2] CRAN (R 4.3.2)
# tzdb                     0.4.0       2023-05-12 [2] CRAN (R 4.3.2)
# utf8                     1.2.4       2023-10-22 [2] CRAN (R 4.3.2)
# vctrs                    0.6.5       2023-12-01 [2] CRAN (R 4.3.2)
# vipor                    0.4.7       2023-12-18 [2] CRAN (R 4.3.2)
# viridis                  0.6.5       2024-01-29 [2] CRAN (R 4.3.2)
# viridisLite              0.4.2       2023-05-02 [2] CRAN (R 4.3.2)

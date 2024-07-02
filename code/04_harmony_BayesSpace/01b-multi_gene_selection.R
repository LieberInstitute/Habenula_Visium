library("spatialLIBD")
# packageVersion("spatialLIBD") # ‘1.15.4’
# library("SpatialExperiment")
# library("scran")
library("tidyverse")
library("here")
library("lobstr")
library("sessioninfo")


# library(dplyr)


source(here("code/04_harmony_BayesSpace", "func_DGE_marker_gene_lists.R")) # Call functions to read paths


### implementing function to filter the data-driven marker gene list from top50r
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
    # cont_colors = viridisLite::turbo(21, direction = 1),
    cont_colors = viridisLite::viridis(21, direction = 1),
    pdf = here(dir_plots, "Marker_genes_exploratory_initial", .y),
    assayname = "logcounts"
))



#######################  Inspect `MHB neuron specific` genes marker #######################

## -----------------------------
# markers.custom$`MHB neuron specific` / `LHB neuron specific` / `mediodorsal thalamus`

suffix_name <- "MHB_neuron_specific.pdf"

# Extract Ensembl ID
lst_Habenula <- rowData(spe)$gene_search[
    rowData(spe)$gene_name %in% markers.custom$`MHB neuron specific`
]

lst_multi_g <- c(
    z_score = paste0("multi_genes_Zs_", suffix_name),
    pca = paste0("multi_genes_PCA_", suffix_name),
    sparsity = paste0("multi_genes_Sp_", suffix_name)
)

print("Ploting multi-gene model for Medial Habenula neuron specific gene markers")

map2(as.vector(names(lst_multi_g)), as.vector(lst_multi_g), ~ vis_grid_gene(
    spe = spe,
    geneid = lst_Habenula,
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



#######################  Inspect `LHB neuron specific` genes marker #######################

## -----------------------------
# markers.custom$`LHB neuron specific` / `mediodorsal thalamus`

suffix_name <- "LHB_neuron_specific.pdf"

# Extract Ensembl ID
lst_Habenula <- rowData(spe)$gene_search[
    rowData(spe)$gene_name %in% markers.custom$`LHB neuron specific`
]

lst_multi_g <- c(
    z_score = paste0("multi_genes_Zs_", suffix_name),
    pca = paste0("multi_genes_PCA_", suffix_name),
    sparsity = paste0("multi_genes_Sp_", suffix_name)
)

print("Ploting multi-gene model for Lateral Habenula neuron specific gene markers")

map2(as.vector(names(lst_multi_g)), as.vector(lst_multi_g), ~ vis_grid_gene(
    spe = spe,
    geneid = lst_Habenula,
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

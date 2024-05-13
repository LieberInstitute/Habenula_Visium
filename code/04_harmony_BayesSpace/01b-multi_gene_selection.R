library("spatialLIBD")
# packageVersion("spatialLIBD") # ‘1.15.4’
# library("SpatialExperiment")
# library("scran")
library("tidyverse")
library("here")
library("lobstr")
library("sessioninfo")


#library(dplyr)


source(here("code/04_harmony_BayesSpace", "func_DGE_marker_gene_lists.R"))       # Call functions to read paths


############################# Load the marker gene list ################################

# We have access to 3 gene markers lists:

# Erik and Top50r putative marker genes merged
markers.custom <- get_erik_and_Hb_markers_genes()          # merged lists
prefix_name <- 'all_gm'                                    # prefix to save matched markers found in the clusters
#markers.custom <- get_bukola_markers_genes_Hb()           # Bukola lists
#prefix_name <- 'erik_gm'  
#markers.custom <- get_Top50r_markers_genes_Hb()           # Top50r lists (putative Hb)
#prefix_name <- 'Top50r_gm'  

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

# set the number of top genes selected for the multi gene analysis
prefix_name <- 'ALL'

####################### Run multi_gene analysis for exploratory purposes ###########################

## Set directory for data and plots

dir_plots <- here("plots", "04_harmony_BayesSpace")
dir_rdata <- here("processed-data", "04_harmony_BayesSpace")


## set path to read RDS object
spe_in_path <- here("processed-data", "04_harmony_BayesSpace", "spe_qc_low_spatialLIBD.rds")


## load Datasets
spe <- readRDS(spe_in_path)
unique(spe$sample_id)
# [1] "V12D07-075_C1" "V13B23-285_A1" "V13B23-285_B1" "V13B23-285_C1"
# [5] "V13B23-285_D1"

cat(" Number of spots:", dim(spe)[2], "\n")

## Set some initials for manage plots
var_height <- 24 # 24/3=8
var_width <- 36 # 36/4=9
var_point_size <- 3.5

# imgData(spe)
# spi <- getImg(spe[18])
# str(spi)
# identical(spi, imgData(spe)$data[[1]])
# plot(imgRaster(spi))


#######################  Inspect WM genes (track)   ####################### 

lst_WM <- c("GFAP", "AQP4", "MBP", "PLP1")

# Extract Ensembl ID
lst_WM <- rowData(spe)$gene_search[
    rowData(spe)$gene_name %in% lst_white_matter_genes
]

## Our list of white matter genes
lstID_WM
# [1] "GFAP; ENSG00000131095" "AQP4; ENSG00000171885" "MBP; ENSG00000197971"
# [4] "PLP1; ENSG00000123560"

## -----------------------------
## in-tissue metrics

lst_multi_g <- c(
  z_score = paste0("multi_genes_Zs_WM_", prefix_name, ".pdf"),
  pca = paste0("multi_genes_PCA_WM_", prefix_name, ".pdf"),
  sparsity = paste0("multi_genes_Sp_WM_", prefix_name, ".pdf")
)

print("Ploting multi-genes for WM gene markers")

map2(as.vector(names(lst_multi_g)), as.vector(lst_multi_g), ~ vis_grid_gene(
  spe = spe,
  geneid = lst_WM,
  multi_gene_method = .x, # z-score, pca, sparcity
  height = var_height, 
  width = var_width,
  point_size = var_point_size,
  #cont_colors = viridisLite::turbo(21, direction = 1),  
  cont_colors = viridisLite::viridis(21, direction = 1),
  #return_plots = TRUE,
  pdf = here(dir_plots, .y),
  assayname = "counts"
))




#######################  Inspect LH genes (track) ####################### 

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
  #cont_colors = viridisLite::turbo(21, direction = 1),  
  cont_colors = viridisLite::viridis(21, direction = 1),
  #return_plots = TRUE,
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
  #cont_colors = viridisLite::turbo(21, direction = 1),  
  cont_colors = viridisLite::viridis(21, direction = 1),
  #return_plots = TRUE,
  pdf = here(dir_plots, .y),
  assayname = "counts"
))


## =============================================================================


modeling_results <- fetch_data(type = "modeling_results")
sce_layer <- fetch_data(type = "sce_layer")
sig_genes <- sig_genes_extract_all(
    n = 400,
    modeling_results = modeling_results,
    sce_layer = sce_layer
)
i_gfap <- subset(sig_genes, gene == "GFAP" &
    test == "WM")$top
i_gfap
set.seed(20200206)
layer_boxplot(
    i = i_gfap,
    sig_genes = sig_genes,
    sce_layer = sce_layer
)



## =============================================================================



## Test white matter (WM) genes

# Define some genes known to be markers for white matter (Tran, Maynard, Spangler, Huuki, Montgomery, Sadashivaiah, Tippani, Barry, Hancock, Hicks, Kleinman, Hyde, Collado-Torres, Jaffe, and Martinowich, 2021). Across five human brain reward circuitry: nucleus accumbens, amygdala, subgenual anterior cingulate cortex, hippocampus, and dorsolateral prefrontal cortex

white_matter_genes <- c("GFAP", "AQP4", "MBP", "PLP1")
white_matter_genes <- rowData(spe)$gene_search[
    rowData(spe)$gene_name %in% white_matter_genes
]

## Our list of white matter genes
white_matter_genes

## WM: Astrocytes and Oligodendrocytes

WM_astro <- c(
    "ENSG00000142611", "ENSG0 0000177133", # PRDM16
    "ENSG00000234377", # OBI1-AS1
    "ENSG00000147509", # RGS20
    "ENSG00000138696", # BMPR1B
    "ENSG00000182902", # SLC25A18
    "ENSG00000164199", # ADGRV1
    "ENSG00000068078", # FGFR3
    "ENSG00000149090", # PAMR1
    "ENSG00000197360", # ZNF98
    "ENSG00000135063"
) # ENTREP1

WM_Oligo <- c(
    "ENSG00000086205", # FOLH1
    "ENSG00000012124", # CD22
    "ENSG00000253877", # LINC01608
    "ENSG00000158865", # SLC5A11
    "ENSG00000105695", # MAG
    # AC012494.1,
    "ENSG00000124920", # MYRF
    "ENSG00000011426", # ANLN
    "ENSG00000204655", # MOG
    "ENSG00000122367"
) # LDB3


WM_genes <- c(WM_astro, WM_Oligo)

# genes = rownames(spe)
# rowData(spe)$gene_search[1:10]

print("Ploting WM pattern in-tissues plot")
sampleID <- unique(spe$sample_id)[3]
# # Subset a specific
# spe_one <- spe[,spe$sample_id == sampleID]
# unique(spe_one$sample_id)

p1 <- vis_gene(
    spe = spe,
    # spe = spe_one,
    sampleid = unique(spe$sample_id)[3],
    geneid = WM_astro, #' MBP',
    spatial = TRUE,
    height = var_height, # 8
    width = var_width, # 9
    point_size = var_point_size,
    return_plots = TRUE,
    # pdf = here(dir_plots, .y),
    assayname = "counts"
)
p1

lobstr::obj_size(spe)

# ## Save object with metrics_qc()
# saveRDS(spe, file.path(dir_rdata, "spe_qc.rds"))


## Save object with metrics_qc()
# saveRDS(spe, file.path(dir_rdata, "spe_qc_low_lib_edge_HighM.rds"))


# ==============================================================================

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

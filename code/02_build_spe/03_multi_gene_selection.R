library("spatialLIBD")
packageVersion("spatialLIBD") # ‘1.15.4’
library("SpatialExperiment")
? # library("scran")
library("tidyverse")
library("here")
library("lobstr")
library("sessioninfo")



# load libraries
library(tidyverse)
library(dplyr)
library(here)

here::here()

# Check if processed_data directory exists, if not create it
if (!dir.exists(here("processed-data/05_DiffExpr_Clustering_CellrangerARC/"))) {
    dir.create(here("processed-data/05_DiffExpr_Clustering_CellrangerARC/"))
}

source(here("code/functions_custom", "remote_DGE_marker_gene_lists.R"))       # Call functions to read paths


#############################           Initials        ################################
############################# Pickup a Marker gene list ################################

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

# set the number of top DGE rows to consider for looking gene markers in the cellranger-arc clusters
n_match_slice <- 20   #10
prefix_name <- paste0(prefix_name, n_match_slice, '.csv')

#############################  Set the DGE list to parse  ################################

## commandArgs scans the arguments which have been supplied when the current R script was invoked (from shell sh)
sample_tmp <- commandArgs(trailingOnly = TRUE)
#sample_tmp <- args[1]
# testing
#sample_tmp <- 'S1_Hb_KDM,human'  # testing HUMAN tissue
#sample_tmp <- 'S2_Hb_KDM,human'  # testing HUMAN tissue
#sample_tmp <- '2_HPC_KDM,human'  # testing HUMAN tissue
#sample_tmp <- '2_HPC_KDM,human'  # testing HUMAN tissue
sample_data = unlist(strsplit(sample_tmp,","))

s_sample <- sample_data[[1]]
s_tissue <- sample_data[[2]]
message('Processing sample: ',s_sample, ' from ', s_tissue, ' tissue.')




## Set directory for data amd plots

dir_plots <- here("plots", "02_build_spe")
dir_rdata <- here("processed-data", "02_build_spe")

## set path to raw and pre-filtered data

spe_in_path <- here("processed-data", "02_build_spe", "spe_qc_low_lib_edge.rds")


## load Datasets

spe <- readRDS(spe_in_path)
unique(spe$sample_id)
# [1] "V12D07-075_C1" "V13B23-285_A1" "V13B23-285_B1" "V13B23-285_C1"
# [5] "V13B23-285_D1"

cat("Initial number of spots:", dim(spe)[2], "\n")

## Set some initials for manage plots
var_height <- 8 # 24/3=8
var_width <- 9 # 36/4=9
var_point_size <- 1.5


## Inspect WM genes (track) tissue

white_matter_genes <- c("GFAP", "AQP4", "MBP", "PLP1")
white_matter_genes <- rowData(spe)$gene_search[
    rowData(spe)$gene_name %in% white_matter_genes
]

## Our list of white matter genes
white_matter_genes
# [1] "GFAP; ENSG00000131095" "AQP4; ENSG00000171885" "MBP; ENSG00000197971"
# [4] "PLP1; ENSG00000123560"

imgData(spe)
spi <- getImg(spe[18])
str(spi)
identical(spi, imgData(spe)$data[[1]])
plot(imgRaster(spi))


## plot 1 gene for 1 sample
vis_gene(
    spe = spe,
    # spe = spe_one,
    sampleid = unique(spe$sample_id)[4],
    geneid = white_matter_genes[1],
    spatial = TRUE,
    height = var_height,
    width = var_width,
    point_size = var_point_size,
    return_plots = TRUE,
    # pdf = here(dir_plots, .y),
    assayname = "counts"
)

## plot n=4 genea for 1 sample
pdf_file <- "in_tissue_multi_genes_WM.pdf"
vis_gene(
    spe = spe,
    sampleid = unique(spe$sample_id)[4],
    geneid = white_matter_genes,
    multi_gene_method = "z_score",
    spatial = TRUE,
    height = var_height,
    width = var_width,
    point_size = var_point_size,
    return_plots = TRUE,
    # pdf = here(dir_plots, pdf_file),
    assayname = "counts"
)

vis_gene(
    spe = spe,
    sampleid = unique(spe$sample_id)[4],
    geneid = white_matter_genes,
    multi_gene_method = "pca",
    spatial = TRUE,
    height = var_height,
    width = var_width,
    point_size = var_point_size,
    return_plots = TRUE,
    # pdf = here(dir_plots, pdf_file),
    assayname = "counts"
)

vis_gene(
    spe = spe,
    sampleid = unique(spe$sample_id)[4],
    geneid = white_matter_genes,
    multi_gene_method = "sparsity",
    spatial = TRUE,
    height = var_height,
    width = var_width,
    point_size = var_point_size,
    return_plots = TRUE,
    # pdf = here(dir_plots, pdf_file),
    assayname = "counts"
)


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

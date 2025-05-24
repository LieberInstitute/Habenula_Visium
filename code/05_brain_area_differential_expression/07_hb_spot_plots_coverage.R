########################################################################
## Make spot plots using the SpatialLIBD app
## Input: spe object with BayesSpace clusters
##
## Authors. CSC
## Data: XXX
## For 60 to 80k spots: $srun --pty --mem=30GB --x11 bash
########################################################################

library("here")
library("spatialLIBD") #[1] ‘1.21.5’
library("purrr")
library("gridExtra")
library("sessioninfo")


## set in/out directories

# spe_dir <- here(
#     "processed-data",
#     "04_harmony_BayesSpace",
#     "spe_harmony_ann.rds"
# )
# re-address path dir
# this spe object is was prepared before to pseudobulk the data to have all the BS meta-data available 
spe_dir <- here(
    "processed-data",
    "05_layer_differential_expression",
    "spe_harmony_ann.rds"
)
dir_plots <- here(
    "plots",
    "05_brain_area_differential_expression",
    "07_hb_spot_plots_coverage"
)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

spe <- readRDS(spe_dir)
spe
colnames(colData(spe))
# [1] "age"        "BayesSpace" "brain_id"   "diagnosis"  "ncells"     "sample_id"  "sex"       
levels(colData(spe)$BayesSpace)
colnames(colData(spe))[grep("BayesSpace_harmony_", colnames(colData(spe)))]

# first make manual verification to identify missing SpD(s) in samples causing visualization issues
table(colData(spe)$sample_id, colData(spe)$BayesSpace)
# 2 samples (V13B23−285_C1, V13B23−285_D1) have a SpD04, but this is missed on the other samples
# Keep only rows where cluster is NOT SpD04 as this SpD (1) have only a few spots, and (2) tat the edge of the tissue in both samples
spe <- spe[, colData(spe)$BayesSpace != "Sp11D04"]
colData(spe)$BayesSpace <- droplevels(colData(spe)$BayesSpace)
levels(colData(spe)$BayesSpace)
table(colData(spe)$sample_id, colData(spe)$BayesSpace)
colData(spe)$BayesSpace_harmony_k11 <- colData(spe)$BayesSpace
# make a spe copy to plot SpatialRegistration data
spe_SR <- spe
len_levels <- length(levels(colData(spe_SR)$BayesSpace))

## =============================================================================

## prepare `spe` with Hb-putative SpD(s) observed in the `Manual Anatomical Annotation` (KDM)

## prepare levels with Hb-putative SpD(s) based on Manual Anatomical Annotation (KDM)
new_labels <- levels(colData(spe_SR)$BayesSpace)
new_labels <- gsub("Sp11D06", "Sp06-HbRNAScope", new_labels)
new_labels <- gsub("Sp11D07", "Sp07-HbRNAScope", new_labels)
new_labels <- gsub("Sp11D11", "Sp11-HbRNAScope", new_labels)
new_labels

# convert to factor with levels in the correct order
colData(spe)$BayesSpace_harmony_k11 <- factor(
    colData(spe)$BayesSpace_harmony_k11,
    levels = levels(colData(spe)$BayesSpace_harmony_k11), 
    labels = new_labels
)
#colData(spe)$BayesSpace_harmony_k11
levels(colData(spe)$BayesSpace_harmony_k11)
# [1] "SpD01"            "SpD02"            "SpD03"            "SpD04"           
# [5] "SpD05"            "Sp06-Putative-Hb" "Sp07-Putative-Hb" "SpD08"           
# [9] "SpD09"            "SpD10"            "Sp11-Putative-Hb"


## =============================================================================

## prepare the `spe_SR` copy with Hb-putative SpD(s) observed in the `SpatialRegistration` Broad - Heatmap

new_labels_SR <- levels(colData(spe_SR)$BayesSpace)
new_labels_SR <- gsub("Sp11D05", "Sp05-HbSpatialR", new_labels_SR)
new_labels_SR <- gsub("Sp11D10", "Sp10-HbSpatialR", new_labels_SR)
new_labels_SR

# convert to factor with levels in the correct order
colData(spe_SR)$BayesSpace_harmony_k11 <- factor(
    colData(spe_SR)$BayesSpace_harmony_k11,
    levels = levels(colData(spe_SR)$BayesSpace_harmony_k11), 
    labels = new_labels_SR
)
levels(colData(spe_SR)$BayesSpace_harmony_k11)
# [1] "Sp11D01"         "Sp11D02"         "Sp11D03"         "Sp05-HbSpatialR"
# [5] "Sp11D06"         "Sp11D07"         "Sp11D08"         "Sp11D09"        
# [9] "Sp10-HbSpatialR" "Sp11D11"
#table(colData(spe_SR)$sample_id, colData(spe_SR)$BayesSpace_harmony_k11)


## =============================================================================

# subset 1 sample by donor for reference - RNAScope Annot
sample_ids_to_keep <- c("V13B23-280_A1", "V13B23-285_B1", "V14F07-340_A1")
# Subset spe to include only the specified sample_ids
spe_subset <- spe[, colData(spe)$sample_id %in% sample_ids_to_keep]
table(colData(spe_subset)$sample_id)
levels(colData(spe_subset)$BayesSpace_harmony_k11)

# subset 1 sample by donor for reference - SpatialRegistration Annots
spe_subset_SR <- spe_SR[, colData(spe_SR)$sample_id %in% sample_ids_to_keep]
table(colData(spe_subset_SR)$sample_id)
levels(colData(spe_subset_SR)$BayesSpace_harmony_k11)

## =============================================================================

## Set some initials for manage spot size in the plots

var_height <- 24 # 24/3=8
var_width <- 26 # 36/4=9
var_point_size <- 1.5

set.seed(07112024)

# lst_order <- sort(unique(spe$sample_id))
# lst_order
# # [1] "V13B23-280_A1" "V13B23-280_B1" "V13B23-280_C1" "V13B23-280_D1"
# # [5] "V13B23-285_A1" "V13B23-285_B1" "V13B23-285_C1" "V13B23-285_D1"
# # [9] "V14F07-340_A1" "V14F07-340_B1" "V14F07-340_C1" "V14F07-340_D1"

## For Kristen talk we pick up BS k=11 given the Hb (manually ann spots) proportion on the samples

## plot SpD(s) with RNAScope Ann
levels(colData(spe)$BayesSpace_harmony_k11)
color_vector <- c("grey", "#b2df8a", "#e41a1c", "darkgreen", "blue",
                  "yellow", "black", "#a65628", "violet", "purple")

p1_lst <- vis_grid_clus(
    spe = spe_subset, #spe,
    clustervar = "BayesSpace_harmony_k11",
    #sample_order = lst_order,
    height = var_height, # 8
    width = var_width, # 9
    point_size = var_point_size,
    sort_clust = FALSE,
    guide_point_size = 5,
    colors = color_vector,
    return_plots = TRUE
)

## plot SpD(s) with SpatialReg Ann
p2_lst <- vis_grid_clus(
    spe = spe_subset_SR,
    clustervar = "BayesSpace_harmony_k11",
    #sample_order = lst_order,
    height = var_height, # 8
    width = var_width, # 9
    point_size = var_point_size,
    sort_clust = FALSE,
    guide_point_size = 5,
    colors = color_vector,
    return_plots = TRUE
)



## plot SpD(s) with HABENULA - RNAScope Ann
color_vector <- c("grey", "grey", "grey", "gray", "blue",
                  "yellow", "gray", "grey", "grey","purple")

p3_lst <- vis_grid_clus(
    spe = spe_subset,
    clustervar = "BayesSpace_harmony_k11",
    #sample_order = lst_order,
    height = var_height,
    width = var_width,
    point_size = var_point_size,
    #pdf = here(dir_plots, paste0("spe_BayesSpace_k11_annotated_Hb_clustergrid.pdf")),
    sort_clust = FALSE,
    guide_point_size = 5,
    colors = color_vector,
    return_plots = TRUE
)

## plot SpD(s) with HABENULA - RNAScope Ann
color_vector <- c("grey", "grey", "grey", "darkgreen", "gray",
                  "gray", "gray", "grey", "violet","gray")
p4_lst <- vis_grid_clus(
    spe = spe_subset_SR,
    clustervar = "BayesSpace_harmony_k11",
    #sample_order = lst_order,
    height = var_height,
    width = var_width,
    point_size = var_point_size,
    #pdf = here(dir_plots, paste0("spe_BayesSpace_k11_annotated_Hb_clustergrid.pdf")),
    sort_clust = FALSE,
    guide_point_size = 5,
    colors = color_vector,
    return_plots = TRUE
)

message("Integrating plots")

all_plots <- c(p1_lst, p2_lst, p3_lst, p4_lst)
length(all_plots)

# Create the PDF
fn <- here(dir_plots, paste0("BayesSpace_k11_Hb_clustergrid.pdf"))
pdf(fn, height = var_height, width = var_width)

grid.arrange(grobs = all_plots, ncol = 3)

dev.off()


message(" Plots DONE!")

# library(slurmjobs)
# slurmjobs::job_single('07_hb_spot_plots_coverage',
#                       create_shell = TRUE, 
#                       memory = '30G',
#                       command = "07_hb_spot_plots_coverage.R",
#                       partition = "katun")


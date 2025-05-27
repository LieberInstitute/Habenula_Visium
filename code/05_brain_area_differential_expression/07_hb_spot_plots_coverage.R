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
library("grid")
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
    "05_brain_area_differential_expression",
    "spe_harmony_ann.rds"
)
dir_plots <- here(
    "plots",
    "05_brain_area_differential_expression",
    "07_hb_spot_plots_coverage"
)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

## load spe with BayesSpace_harmony data
spe <- readRDS(spe_dir)
spe

# Verification
colnames(colData(spe))
colnames(colData(spe))[grep("BayesSpace_harmony_", colnames(colData(spe)))]
levels(colData(spe)$BayesSpace)

## set levels to desired 'BayesSpace_harmony_k', and format as desired 
# BayesSpace clustering 
k=11
k_nice <- sprintf("%02d", k)

message("Processing spotPlots for BS k=", k_nice)

# Set the levels of colData(spe)$BayesSpace to be the levels of BayesSpace_harmony_k11
colData(spe)$BayesSpace <- colData(spe)$BayesSpace_harmony_k11
# check levels of both
unique(colData(spe)$BayesSpace_harmony_k11)
unique(colData(spe)$BayesSpace)
levels(colData(spe)$BayesSpace)
# NULL

## format level as desired
colData(spe)$BayesSpace <- paste0("Sp", k_nice, 
                                  "D", sprintf("%02d", as.numeric(as.character(colData(spe)$BayesSpace)))
                                  )
unique(colData(spe)$BayesSpace)
colData(spe)$BayesSpace <- factor(colData(spe)$BayesSpace)
levels(colData(spe)$BayesSpace)
# [1] "Sp11D01" "Sp11D02" "Sp11D03" "Sp11D04" "Sp11D05" "Sp11D06" "Sp11D07"
# [8] "Sp11D08" "Sp11D09" "Sp11D10" "Sp11D11"

## make manual verification to identify missing SpD(s) in samples causing visualization issues
spots_by_domain_table <- table(colData(spe)$sample_id, colData(spe)$BayesSpace)
spots_by_domain_table
# Get the table of counts for each level in colData(spe)$BayesSpace
bayes_space_counts <- colSums(spots_by_domain_table)
bayes_space_counts
# Identify levels with less than 20 counts
levels_to_remove <- names(bayes_space_counts[bayes_space_counts < 20])
message("Spatial Domains with less than 20 counts: ", paste(levels_to_remove, collapse = ", "))

# Remove those levels from BayesSpace by excluding them
colData(spe)$BayesSpace <- factor(colData(spe)$BayesSpace, 
                                  levels = setdiff(levels(colData(spe)$BayesSpace), levels_to_remove))
# Check the new levels of BayesSpace
levels(colData(spe)$BayesSpace)
# [1] "Sp11D01" "Sp11D02" "Sp11D03" "Sp11D05" "Sp11D06" "Sp11D07" "Sp11D08"
# [8] "Sp11D09" "Sp11D10" "Sp11D11"

## avoid manual evaluation
# # 2 samples (V13B23−285_C1, V13B23−285_D1) have a SpD04, but this is missed on the other samples
# # Keep only rows where cluster is NOT SpD04 as this SpD (1) have only a few spots, and (2) tat the edge of the tissue in both samples
# spe <- spe[, colData(spe)$BayesSpace != "Sp11D04"]
# colData(spe)$BayesSpace <- droplevels(colData(spe)$BayesSpace)
# levels(colData(spe)$BayesSpace)
# table(colData(spe)$sample_id, colData(spe)$BayesSpace)
# colData(spe)$BayesSpace_harmony_k11 <- colData(spe)$BayesSpace

# make a spe copy to plot SpatialRegistration data vs Hb RNScope taxonomy (KDM)
spe_SR <- spe
spe_SR
len_levels <- length(levels(colData(spe_SR)$BayesSpace)) # sould be NULL to add the SR SpD
message("Spatial Domains after removing SpD(s) with less than 20 counts: ", len_levels)

## =============================================================================
## subset 1 sample by donor for reference

sample_ids_to_keep <- c("V13B23-280_A1", "V13B23-285_B1", "V14F07-340_A1")
# Subset spe to include only the specified sample_ids
spe_subset <- spe[, colData(spe)$sample_id %in% sample_ids_to_keep]
table(colData(spe_subset)$sample_id)
levels(colData(spe_subset)$BayesSpace)

# subset 1 sample by donor for reference - SpatialRegistration Annots
spe_subset_SR <- spe_SR[, colData(spe_SR)$sample_id %in% sample_ids_to_keep]
table(colData(spe_subset_SR)$sample_id)
levels(colData(spe_subset_SR)$BayesSpace)

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

## Plot ALL domains for reference

## plot SpD(s) with RNAScope Ann
levels(colData(spe)$BayesSpace)
color_vector <- c("grey", "#b2df8a", "#e41a1c", "darkgreen", "blue",
                  "yellow", "black", "#a65628", "violet", "gold")

p1_lst <- vis_grid_clus(
    spe = spe_subset, #spe,
    clustervar = "BayesSpace",
    #sample_order = lst_order,
    height = var_height, # 8
    width = var_width, # 9
    point_size = var_point_size,
    sort_clust = FALSE,
    guide_point_size = 5,
    colors = color_vector,
    return_plots = TRUE
)

# ## plot SpD(s) with SpatialReg Ann
# p2_lst <- vis_grid_clus(
#     spe = spe_subset_SR,
#     clustervar = "BayesSpace",
#     #sample_order = lst_order,
#     height = var_height, # 8
#     width = var_width, # 9
#     point_size = var_point_size,
#     sort_clust = FALSE,
#     guide_point_size = 5,
#     colors = color_vector,
#     return_plots = TRUE
# )


## =============================================================================

## prepare `spe` with Hb-putative SpD(s) observed in the `Manual Anatomical Annotation` (KDM)

## prepare levels with Hb-putative SpD(s) based on Manual Anatomical Annotation (KDM)
new_labels <- levels(colData(spe_subset)$BayesSpace)
new_labels <- gsub("Sp11D06", "Sp06-HbRNAScope", new_labels)
new_labels <- gsub("Sp11D07", "Sp07-HbRNAScope", new_labels)
new_labels <- gsub("Sp11D11", "Sp11-HbRNAScope", new_labels)
new_labels

## convert to factor with levels in the correct order
colData(spe_subset)$BayesSpace <- factor(
    colData(spe_subset)$BayesSpace,
    levels = levels(colData(spe_subset)$BayesSpace),
    labels = new_labels
)
levels(colData(spe_subset)$BayesSpace)
# [1] "Sp11D01"         "Sp11D02"         "Sp11D03"         "Sp11D05"        
# [5] "Sp06-HbRNAScope" "Sp07-HbRNAScope" "Sp11D08"         "Sp11D09"        
# [9] "Sp11D10"         "Sp11-HbRNAScope"

## =============================================================================

## prepare the `spe_SR` copy with Hb-putative SpD(s) observed in the `SpatialRegistration` Broad - Heatmap

new_labels_SR <- levels(colData(spe_subset_SR)$BayesSpace)
new_labels_SR <- gsub("Sp11D06", "Sp06-HbSpatialR", new_labels_SR)
new_labels_SR <- gsub("Sp11D11", "Sp11-HbSpatialR", new_labels_SR)
new_labels_SR

## convert to factor with levels in the correct order
colData(spe_subset_SR)$BayesSpace <- factor(
    colData(spe_subset_SR)$BayesSpace,
    levels = levels(colData(spe_subset_SR)$BayesSpace),
    labels = new_labels_SR
)
levels(colData(spe_subset_SR)$BayesSpace)

## =============================================================================

## plot SpD(s) with HABENULA - RNAScope Ann
color_vector <- c("grey", "grey", "grey", "gray", "blue",
                  "yellow", "gray", "grey", "grey","gold")

p3_lst <- vis_grid_clus(
    spe = spe_subset,
    clustervar = "BayesSpace",
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

## plot SpD(s) with HABENULA - SpatialRegistration Ann
color_vector <- c("grey", "grey", "grey", "gray", "blue",
                  "gray", "gray", "grey", "gray","gold")
p4_lst <- vis_grid_clus(
    spe = spe_subset_SR,
    clustervar = "BayesSpace",
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

# Create the PDF
fn <- here(dir_plots, paste0("BayesSpace_k", k_nice, "_Hb_clustergrid.pdf"))

# Combine the lists and add titles to each group
all_plots <- grid.arrange(
    gridExtra::arrangeGrob(grobs = p1_lst, ncol = 3, top = textGrob("All SpD(s) for one sample by donor", 
                                                                    gp = gpar(fontsize = 24, fontface = "bold"))),
    gridExtra::arrangeGrob(grobs = p3_lst, ncol = 3, top = textGrob("SpD(s) identified: Hb Annatomical RNAScope Annotations", 
                                                                    gp = gpar(fontsize = 24, fontface = "bold"))),
    gridExtra::arrangeGrob(grobs = p4_lst, ncol = 3, top = textGrob("SpD(s) identified: SpatialRegistration Correlations", 
                                                                    gp = gpar(fontsize = 24, fontface = "bold"))),
    nrow = 3,
    top = textGrob(paste("Overall SpD(s) for k=", k_nice, "\n"), gp = gpar(fontsize = 28, fontface = "bold"))
)

# all_plots <- c(p1_lst, p3_lst, p4_lst)
# length(all_plots)

pdf(fn, height = var_height, width = var_width)

grid.draw(all_plots)

dev.off()


message(" Plots DONE!")

# library(slurmjobs)
# slurmjobs::job_single('07_hb_spot_plots_coverage',
#                       create_shell = TRUE, 
#                       memory = '30G',
#                       command = "07_hb_spot_plots_coverage.R",
#                       partition = "katun")


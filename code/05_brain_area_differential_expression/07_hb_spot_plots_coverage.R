########################################################################
## Make spot plots using the SpatialLIBD app
## Input: spe object with BayesSpace clusters
## Output: some stats in logs and spotPlot with both 
##.        - (A) SpD(s) selected based on Hb-Anatomical annotations
##.        - (B) SpD(s) selected based on Spatial-Registration annotations
##
## Authors. CSC
## Data: XXX
## For 60 to 80k spots: $srun --pty --mem=30GB --x11 bash
########################################################################

library("here")
library("spatialLIBD") #[1] ‘1.21.5’
library("purrr")
library("grid")
library("viridis")
library("RColorBrewer")
library("gridExtra")
library("sessioninfo")

args = commandArgs(trailingOnly = TRUE)
k <- as.integer(args[2])

if (is.na(k)) {
    k <- 11
}

## set in/out directories
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
k_nice <- sprintf("%02d", k)
BayesSpace_harmony_k <- paste0("BayesSpace_harmony_k", k_nice)

message("Processing spotPlots for BS k=", k_nice)

# Set the levels of colData(spe)$BayesSpace to be the levels of BayesSpace_harmony_k
colData(spe)$BayesSpace <- colData(spe)[[BayesSpace_harmony_k]]
# check levels of both
unique(colData(spe)[[BayesSpace_harmony_k]])
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


## =============================================================================
## Additional analysis
## Get count and percentages for BS k of interest

if (k==11 || k==15 || k==20 || k==28) {
    
    if (k==11) { domains_of_interest <- c("Sp11D06", "Sp11D11") }
    if (k==15) { domains_of_interest <- c("Sp15D06", "Sp15D10", "Sp15D14") }
    if (k==20) { domains_of_interest <- c("Sp20D06", "Sp20D08", "Sp20D16", "Sp20D20") }
    if (k==28) { domains_of_interest <- c("Sp28D05", "Sp28D10", "Sp28D11", "Sp28D20", "Sp28D27") }
    
    message("Hb Spatial Domians in k", k_nice)
    #domains_of_interest <- c("Sp28D05", "Sp28D10", "Sp28D11", "Sp28D20", "Sp28D27")
    bayes_space_counts <- colSums(spots_by_domain_table)[domains_of_interest] 
    percentage <- round((bayes_space_counts * 100) / sum(spots_by_domain_table), 2)
    result_table <- data.frame(
        Domain = domains_of_interest,
        Count = bayes_space_counts,
        Percentage = paste(percentage, "%")
    )
    rownames(result_table) <- NULL
    # Calculate the total percentage and add it as a row at the end of the table
    total_percentage <- round(sum(percentage), 2)
    result_table <- rbind(result_table, data.frame(Domain = "Total", Count = sum(bayes_space_counts), Percentage = paste(total_percentage, "%")))
    message("Counts and percents for BayesSpace k", k_nice)
    print(result_table, row.names = FALSE)
    # Domain Count Percentage
    # Sp28D05   501      1.5 %
    # Sp28D10   519     1.55 %
    # Sp28D11   472     1.41 %
    # Sp28D20   410     1.23 %
    # Sp28D27   350     1.05 %
    # Total  2252     6.74 %
}

## =============================================================================

# Identify levels with less than 20 counts
levels_to_remove <- names(bayes_space_counts[bayes_space_counts < 20])

# Remove those levels from BayesSpace by excluding them
if (length(levels_to_remove) > 0) {
    message("Spatial Domains with less than 20 counts: ", paste(levels_to_remove, collapse = ", "))
    colData(spe)$BayesSpace <- factor(colData(spe)$BayesSpace, 
                                      levels = setdiff(levels(colData(spe)$BayesSpace), levels_to_remove))
}
# Check the new levels of BayesSpace
levels(colData(spe)$BayesSpace)
# [1] "Sp11D01" "Sp11D02" "Sp11D03" "Sp11D05" "Sp11D06" "Sp11D07" "Sp11D08"
# [8] "Sp11D09" "Sp11D10" "Sp11D11"

# make a spe copy to plot SpatialRegistration data vs Hb RNScope taxonomy (KDM)
spe_SR <- spe
spe_SR
len_levels <- length(levels(colData(spe_SR)$BayesSpace)) # should be NULL to add the SR SpD
message("Spatial Domains after removing SpD(s) with less than 20 counts: ", len_levels)

## =============================================================================
## subset 1 sample by donor for reference

sample_ids_to_keep <- c("V13B23-280_A1", "V13B23-285_B1", "V14F07-340_A1")
# Subset spe to include only the specified sample_ids
spe_subset <- spe[, colData(spe)$sample_id %in% sample_ids_to_keep]
table(colData(spe_subset)$sample_id)
levels(colData(spe_subset)$BayesSpace)

# subset 1 sample by donor for reference - SpatialRegistration Annots
spe_subset_SR <- spe_subset
# spe_subset_SR <- spe_SR[, colData(spe_SR)$sample_id %in% sample_ids_to_keep]
# table(colData(spe_subset_SR)$sample_id)
# levels(colData(spe_subset_SR)$BayesSpace)

## =============================================================================

## Make color vector for BS k of interest

## Define the function to generate the color lists with additional custom position inputs
# - creates 4 vectors of colors for handle custom palettes on the spot Plots
# - Includes: colors for ALL the clusters, for the RNAScope custom clusters and for the SpatialReg custom clusters

generate_color_lists <- function(num_colors, custom_RNAScope_SpD, custom_SpatialReg_SpD) {
    #library(RColorBrewer)
    # Generate a vector of colors from the viridis palette based on the input number
    #random_colors <- viridis(num_colors)
    random_colors <- c(
        "red", "blue", "green", "yellow", "orange", 
        "purple", "brown", "pink", "cyan", "magenta", 
        "black", "#24FF24", "darkblue", "darkgreen", "lightblue", 
        "darkorange", "violet", "gold", "darkred", "indianred"
    )
    if (num_colors>20) {
        # Additional 20 colors to be added to the vector for k(s)>20
        additional_colors <- c(
            "lightgreen", "darkgray", "lightgray", "#290AD8", "#1E8E99", 
            "maroon", "midnightblue", "darkviolet", "steelblue", "slateblue", 
            "#A50021", "#FFFF6D", "peachpuff", "plum", "yellowgreen", 
            "turquoise", "chocolate", "firebrick", "salmon", "#490092"
        )
        random_colors <- c(random_colors, additional_colors)
    }
    
    random_colors_SR <- random_colors
    
    # Create gray vector for RNAScope
    grey_colors_RNAScope <- rep("grey", num_colors)
    grey_colors_SpatialReg <- grey_colors_RNAScope
    
    # Adjust grey vector colors based on custom_RNAScope_SpD input (position adjustment by -1)
    grey_colors_RNAScope[custom_RNAScope_SpD] <- random_colors[custom_RNAScope_SpD]
    
    # Adjust grey vector for SpatialReg based on custom_SpatialReg_SpD input (position adjustment by -1)
    grey_colors_SpatialReg[custom_SpatialReg_SpD] <- random_colors_SR[custom_SpatialReg_SpD]
    
    # Return the four lists
    return(list(
        random_colors = random_colors,
        #random_colors_SR = random_colors_SR,
        grey_colors_RNAScope = grey_colors_RNAScope,
        grey_colors_SpatialReg = grey_colors_SpatialReg
    ))
}


# Note: SpD(s) removed break the continuous of the color palette driving to color issues when rendering the spotPlot visualization
# depending on the location of levels removed, these are handled manually 
# g.e. if SpD5 is removed, I move back the color vector for grey(s) after the SpD5 position
    
if (k_nice=="03") {
    # for handle palette of colors
    custom_RNAScope_SpD <- 3
    custom_SpatialReg_SpD <- 3
    # for handle legends labels
    RNAScope_SpD <- 3
    SpatialReg_SpD <- 3
    
} else if (k_nice=="11") {
    # for handle palette of colors
    custom_RNAScope_SpD <- c(6, 7, 11) - length(levels_to_remove) 
    custom_SpatialReg_SpD <- c(6, 11) - length(levels_to_remove)
    # for handle legends labels
    RNAScope_SpD <- c(6, 7, 11)
    SpatialReg_SpD <- c(6, 11)
    
} else if (k_nice=="15") {
    # for handle palette of colors
    custom_RNAScope_SpD <- c(6, 8, 10, 14) - length(levels_to_remove)
    custom_SpatialReg_SpD <- c(6, 10, 14) - length(levels_to_remove)
    # for handle legends labels
    RNAScope_SpD <- c(6, 8, 10, 14)
    SpatialReg_SpD <- c(6, 10, 14)
    
} else if (k_nice=="20") {
    # for handle palette of colors
    custom_RNAScope_SpD <- c(6, 7, 8, 15, 16, 17, 19) - length(levels_to_remove)
    custom_SpatialReg_SpD <- c(6, 8, 16, 19)  - length(levels_to_remove)
    # for handle legends labels
    RNAScope_SpD <- c(6, 7, 8, 15, 16, 17, 19)
    SpatialReg_SpD <- c(6, 8, 16, 19)
    
} else if (k_nice=="28") {
    #levels_to_remove # [1] "Sp28D18" "Sp28D26"
    # for handle palette of colors
    custom_RNAScope_SpD <- c(5, 10, 11, 16, c(20) - 1, c(27, 28) - 2) # - length(levels_to_remove); this do not apply because of SpD removed aren't continuous
    custom_SpatialReg_SpD <- c(5, 10, 11, c(20) - 1, c(27) - 2)
    # for handle legends labels
    RNAScope_SpD <- c(5, 10, 11, 16, 20, 27, 28)
    SpatialReg_SpD <- c(5, 10, 11, 20, 27)
    
}    


# Call the function with the new pallete of colors
n_colors <- length(levels(colData(spe_subset)$BayesSpace))
color_lists <- generate_color_lists(n_colors, custom_RNAScope_SpD, custom_SpatialReg_SpD)

# verification
print(color_lists$random_colors)
#print(color_lists$random_colors_SR)
print(color_lists$grey_colors_RNAScope)
print(color_lists$grey_colors_SpatialReg)


## =============================================================================

## Set some initials for manage spot size in the plots

var_height <- 24 # 24/3=8
var_width <- 26 # 36/4=9
var_point_size <- 2
var_guide_point_size = 5
if (k > 15) { var_guide_point_size = 3 }

set.seed(07112024)

# lst_order <- sort(unique(spe$sample_id))
# lst_order
# # [1] "V13B23-280_A1" "V13B23-280_B1" "V13B23-280_C1" "V13B23-280_D1"
# # [5] "V13B23-285_A1" "V13B23-285_B1" "V13B23-285_C1" "V13B23-285_D1"
# # [9] "V14F07-340_A1" "V14F07-340_B1" "V14F07-340_C1" "V14F07-340_D1"

## Plot ALL domains for reference

## plot SpD(s) with RNAScope Ann
levels(colData(spe)$BayesSpace)
## set color_vector as corresponds 
color_vector <- color_lists$random_colors

p1_lst <- vis_grid_clus(
    spe = spe_subset,
    clustervar = "BayesSpace",
    #sample_order = lst_order,
    height = var_height, # 8
    width = var_width, # 9
    point_size = var_point_size,
    sort_clust = FALSE,
    guide_point_size = var_guide_point_size,
    colors = color_lists$random_colors,
    return_plots = TRUE
)


## =============================================================================

## prepare `spe` with Hb-putative SpD(s) observed in the `Manual Anatomical Annotation` (KDM)

new_labels <- levels(colData(spe_subset)$BayesSpace)
new_labels_SR <- new_labels
RNAScope_SpD <- sprintf("%02d", RNAScope_SpD)
old_labels <- unlist(map(RNAScope_SpD, ~ paste0("Sp", k_nice, "D", .x)))

for (i in seq_along(RNAScope_SpD)) {
    #new_labels <- gsub(old_labels[i], paste0(old_labels[i], "-HbRNAScope"), new_labels)
    # Shorter legend labels for visualization purposes. Extract the substring starting from the second "D"
    substring_after_second_D <- sub(".*D(.*)", "D\\1", old_labels[i])
    new_labels <- gsub(old_labels[i], paste0(substring_after_second_D, "-HbRNAScope"), new_labels)
}
     
new_labels

## convert to factor with levels in the correct order
colData(spe_subset)$BayesSpace <- factor(
    colData(spe_subset)$BayesSpace,
    levels = levels(colData(spe_subset)$BayesSpace),
    labels = new_labels
)
levels(colData(spe_subset)$BayesSpace)


## =============================================================================

## prepare the `spe_SR` copy with Hb-putative SpD(s) observed in the `SpatialRegistration` Broad - Heatmap
   
new_labels_SR #<- levels(colData(spe)$BayesSpace)
SpatialReg_SpD <- sprintf("%02d", SpatialReg_SpD)
old_labels <- unlist(map(SpatialReg_SpD, ~ paste0("Sp", k_nice, "D", .x)))

for (i in seq_along(SpatialReg_SpD)) {
    #new_labels_SR <- gsub(old_labels[i], paste0(old_labels[i], "-HbSpatialR"), new_labels_SR)
    # Shorter legend labels for visualization purposes. Extract the substring starting from the second "D"
    substring_after_second_D <- sub(".*D(.*)", "D\\1", old_labels[i])
    new_labels_SR <- gsub(old_labels[i], paste0(substring_after_second_D, "-HbSpatialR"), new_labels_SR)
}
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

## set color_vector as correspond 
color_vector <- color_lists$grey_colors_RNAScope

p3_lst <- vis_grid_clus(
    spe = spe_subset,
    clustervar = "BayesSpace",
    #sample_order = lst_order,
    height = var_height,
    width = var_width,
    point_size = var_point_size,
    #pdf = here(dir_plots, paste0("spe_BayesSpace_k11_annotated_Hb_clustergrid.pdf")),
    sort_clust = FALSE,
    guide_point_size = var_guide_point_size,
    colors = color_lists$grey_colors_RNAScope,
    return_plots = TRUE
)

## plot SpD(s) with HABENULA - SpatialRegistration Ann

## set color_vector as correspond 
color_vector <- color_lists$grey_colors_SpatialReg

p4_lst <- vis_grid_clus(
    spe = spe_subset_SR,
    clustervar = "BayesSpace",
    #sample_order = lst_order,
    height = var_height,
    width = var_width,
    point_size = var_point_size,
    #pdf = here(dir_plots, paste0("spe_BayesSpace_k11_annotated_Hb_clustergrid.pdf")),
    sort_clust = FALSE,
    guide_point_size = var_guide_point_size,
    colors = color_lists$grey_colors_SpatialReg,
    return_plots = TRUE
)


message("Integrating plots")

# Create the PDF
fn <- here(dir_plots, paste0("BayesSpace_k", k_nice, "_Hb_clustergrid_RNAScope_vs_AnatomicalAnn.pdf"))

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
#      
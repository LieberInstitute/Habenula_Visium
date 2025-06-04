########################################################################
## Compare two BayesSpace k(s)
## Input: spe object with BayesSpace clusters
## Output: spot plot with cluster agreement btw k20 & k28 
##
## Authors. CSC
## Data: June 03, 2025
## For 60 to 80k spots: $srun --pty --mem=15GB --x11 bash
########################################################################

library("here")
library("spatialLIBD")
library("purrr")
library("ggplot2")
library("patchwork")
library("dplyr")
library("sessioninfo")


## set in/out directories
spe_dir <- here(
    "processed-data",
    "05_brain_area_differential_expression",
    "spe_harmony_ann.rds"
)
dir_plots <- here(
    "plots",
    "05_brain_area_differential_expression",
    "09_comparing_BayesSpace_clusters"
)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

## load spe with BayesSpace_harmony data
spe <- readRDS(spe_dir)
spe

# Verification
colnames(colData(spe))
colnames(colData(spe))[grep("BayesSpace_harmony_", colnames(colData(spe)))]


## Add cluster agreement category: Only_k20, Only_k28, Both, or None

# define Habenula-relevant SpD(s) for k=20 and k=28
habenula_k20 <- c(6, 8, 16, 19)
habenula_k28 <- c(5, 10, 11, 20, 27)

## filter specific SpD(s) for each k and add categories
idx_habenula <- spe$BayesSpace_harmony_k20 %in% habenula_k20 |
    spe$BayesSpace_harmony_k28 %in% habenula_k28
length(idx_habenula)
## check overlaps and counts for each k
table(idx_habenula)
table(spe$BayesSpace_harmony_k20 %in% habenula_k20)
table(spe$BayesSpace_harmony_k28 %in% habenula_k28)

## Pre-fill all values with "None"
colData(spe)$cluster_agreement <- factor("None", levels = c("Only_k20", "Only_k28", "Both", "None"))
levels(colData(spe)$cluster_agreement)
table(spe$cluster_agreement)
# Only_k20 Only_k28     Both     None 
#       0        0        0    33409

## Subset values using the same index for both
k20 <- as.character(spe$BayesSpace_harmony_k20[idx_habenula])
k28 <- as.character(spe$BayesSpace_harmony_k28[idx_habenula])

## make categories
## define agreement logic
labels <- case_when(
    k20 %in% habenula_k20 & k28 %in% habenula_k28 ~ "Both",
    k20 %in% habenula_k20 & !(k28 %in% habenula_k28) ~ "Only_k20",
    !(k20 %in% habenula_k20) & k28 %in% habenula_k28 ~ "Only_k28",
    TRUE ~ "None"
)

# Assign labels back only to matching rows
colData(spe)$cluster_agreement[idx_habenula] <- factor(labels, levels = c("Only_k20", "Only_k28", "Both", "None"))
table(colData(spe)$cluster_agreement)
# Only_k20 Only_k28     Both     None 
#       81       39     2213    31076

# ## Subset spe to selected samples for visualization ============================
# 
# sample_ids_to_keep <- c("V13B23-280_A1", "V13B23-285_B1", "V14F07-340_A1")
# 
# # Subset the SpatialExperiment object
# spe <- spe[, colData(spe)$sample_id %in% sample_ids_to_keep]
# # Confirm sample IDs and levels
# table(colData(spe)$sample_id)
# levels(colData(spe)$cluster_agreement)
# ##  ============================ ============================ ==================

## prepare data for plotting

# Extract spatial coordinates and combine with colData
df_plot <- cbind(
    as.data.frame(colData(spe)),
    spatialCoords(spe)
)

## spatial plot of cluster agreement

# ggplot(df_plot, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres)) +
#     geom_point(aes(color = cluster_agreement), size = 0.5) +
#     scale_color_manual(
#         values = c(
#             "Only_k20" = "#1b9e77",
#             "Only_k28" = "#d95f02",
#             "Both"     = "#7570b3",
#             "None"     = "#cccccc"
#         ),
#         name = "Cluster Agreement"
#     ) +
#     scale_y_reverse() +  # Match spatial layout
#     coord_fixed() +
#     facet_wrap(~ sample_id) +
#     theme_minimal() +
#     labs(
#         title = "Spatial Agreement of BayesSpace k=20 and k=28 Clusters",
#         x = "Pixel Column",
#         y = "Pixel Row"
#     )



var_height <- 24 # 24/3=8
var_width <- 26 # 36/4=9
var_point_size <- 2
var_guide_point_size = 3

set.seed(07112024)

p1_lst <- vis_grid_clus(
    spe = spe,
    clustervar = "cluster_agreement",
    height = var_height,
    width = var_width,
    point_size = var_point_size,
    sort_clust = FALSE,
    guide_point_size = var_guide_point_size,
    colors = c("red", "gold", "darkgreen", "grey"), ## Only_k20 Only_k28     Both     None 
    return_plots = TRUE
)

p1_integrated <- wrap_plots(p1_lst, nrow = 3, ncol = 4)

f_name <- here(dir_plots, "BS_k20_k28_comparison.pdf")
pdf(f_name, height = var_height, width = var_width)

print(p1_integrated)

dev.off()



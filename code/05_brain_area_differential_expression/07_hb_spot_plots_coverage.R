library("here")
library("spatialLIBD") #[1] ‘1.21.4’
library("purrr")
library("sessioninfo")


## set in/out directories

spe_dir <- here(
    "processed-data",
    "04_harmony_BayesSpace",
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

# current SpDs annotated in pseudobulk spe
# k13
colnames(colData(spe))[grep("BayesSpace_harmony_", colnames(colData(spe)))]
levels(colData(spe)$BayesSpace)
# [1] "Sp13D01.Oligo"      "Sp13D02.Endo"       "Sp13D03.Endo"
# [4] "Sp13D04.Astrocyte"  "Sp13D05.Inhib.Thal" "Sp13D06.Astrocyte"
# [7] "Sp13D07.Astrocyte"  "Sp13D08.OPC_Astroc" "Sp13D09.Astrocyte_"
# [10] "Sp13D10.Oligo"      "Sp13D11.Habenula"   "Sp13D12.Oligo"
# [13] "Sp13D13.Oligo"

# prepare new levels to plot by other BS k of interest
levels(colData(spe)$BayesSpace_harmony_k11)
colData(spe)$BayesSpace_harmony_k11 <- factor(colData(spe)$BayesSpace_harmony_k11)
levels(colData(spe)$BayesSpace_harmony_k11)

# Set custom levels
new_labels <- paste0("SpD", sprintf("%02d", seq(1:11)))
new_labels <- gsub("SpD06", "Sp06-Putative-Hb", new_labels)
new_labels <- gsub("SpD07", "Sp07-Putative-Hb", new_labels)
new_labels <- gsub("SpD11", "Sp11-Putative-Hb", new_labels)
new_labels

# convert to factor with levels in the correct order
colData(spe)$BayesSpace_harmony_k11 <- factor(
    colData(spe)$BayesSpace_harmony_k11,
    levels = 1:11,
    labels = new_labels
)
levels(colData(spe)$BayesSpace_harmony_k11)
# [1] "SpD01"            "SpD02"            "SpD03"            "SpD04"           
# [5] "SpD05"            "Sp06-Putative-Hb" "Sp07-Putative-Hb" "SpD08"           
# [9] "SpD09"            "SpD10"            "Sp11-Putative-Hb"


## Set some initials for manage spot size in the plots

var_height <- 24 # 24/3=8
var_width <- 26 # 36/4=9
var_point_size <- 1.5

set.seed(07112024)

lst_order <- sort(unique(spe$sample_id))
# [1] "V13B23-280_A1" "V13B23-280_B1" "V13B23-280_C1" "V13B23-280_D1"
# [5] "V13B23-285_A1" "V13B23-285_B1" "V13B23-285_C1" "V13B23-285_D1"
# [9] "V14F07-340_A1" "V14F07-340_B1" "V14F07-340_C1" "V14F07-340_D1"

# vis_grid_clus(
#     spe = spe,
#     clustervar = "BayesSpace",
#     sample_order = lst_order,
#     height = var_height, # 8
#     width = var_width, # 9
#     point_size = var_point_size,
#     pdf = here(dir_plots, paste0("spe_BayesSpace_k", k ,"_annotated_clustergrid.pdf")),
#     sort_clust = FALSE,
#     guide_point_size = 5,
#     colors = c("grey", "#b2df8a", "#e41a1c", "#377eb8", "#4daf4a",
#                "#ff7f00", "black", "#a65628", "#999999", "blue",
#                "purple", "gold")
# )

## For Kristen talk we pick up BS k=11 given the Hb (manually ann spots) proportion on the samples
levels(colData(spe)$BayesSpace_harmony_k11)
# [1] "SpD01"            "SpD02"            "SpD03"            "SpD04"           
# [5] "SpD05"            "Sp06-Putative-Hb" "Sp07-Putative-Hb" "SpD08"           
# [9] "SpD09"            "SpD10"            "Sp11-Putative-Hb"

## test: remove color element for SpD04
color_vector <- c("grey", "grey", "grey", "grey", "blue", "yellow", "grey", "grey", "grey","purple")

vis_grid_clus(
    spe = spe,
    clustervar = "BayesSpace_harmony_k11",
    sample_order = lst_order,
    height = var_height,
    width = var_width,
    point_size = var_point_size,
    pdf = here(dir_plots, paste0("spe_BayesSpace_k11_annotated_Hb_clustergrid.pdf")),
    sort_clust = FALSE,
    guide_point_size = 5,
    colors = color_vector
)


# # Trying to add a fake domain for visualization purposes on he spot plots
# samples <- unique(colData(spe)$sample_id)
# 
# for (s in unique(colData(spe)$sample_id)) {
#     spe_sub <- spe[, colData(spe)$sample_id == s]
#     
#     # Ensure correct levels in factor
#     colData(spe_sub)$BayesSpace_harmony_k11 <- factor(
#         colData(spe_sub)$BayesSpace_harmony_k11,
#         levels = 1:11
#     )
#     
#     # Identify missing levels
#     missing_lvls <- setdiff(1:11, unique(as.integer(colData(spe_sub)$BayesSpace_harmony_k11)))
#     
#     if (length(missing_lvls) > 0) {
#         # Use the first spot to clone structure
#         dummy_template <- spe_sub[, 1]
#         
#         # Create dummy columns, one for each missing level
#         for (lvl in missing_lvls) {
#             dummy_col <- dummy_template
#             
#             # Replace metadata with dummy values
#             colData(dummy_col)$BayesSpace_harmony_k11 <- factor(lvl, levels = 1:11)
#             
#             # Optional: shift spatial coordinates far off the plot area
#             spatialCoords(dummy_col)[,1] <- -9999
#             spatialCoords(dummy_col)[,2] <- -9999
#             
#             # Combine dummy column
#             spe_sub <- cbind(spe_sub, dummy_col)
#         }
#     }
#     
#     # Optional: plot using BayesSpace
#     # p <- spatialPlot(spe_sub, label = "BayesSpace_harmony_k11")
#     # plots[[s]] <- p
# }
# 
# color_vector <- c("grey", "grey", "grey", "grey", "grey", "blue", "yellow", "grey", "grey", "grey","purple")


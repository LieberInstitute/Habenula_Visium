
## load spe data

spe_in <- here("processed-data", "04_harmony_BayesSpace", "spe_harmony_ann.rds")
spe <- readRDS(spe_in)

colnames(colData(spe))
clustervar <- colData(spe)$BayesSpace
dir_plots <- here("plots","05_brain_area_differential_expression", "07_hb_spot_plots_coverage")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

unique(spe$sample_id)

## Subset to two samples of interest and obtain the plot list
## Get number of TRUE spots in tissue
in_tissue_spots <- map(unique(spe$sample_id), ~ summary(spe$in_tissue[spe$sample_id == .x] == TRUE))


## Set some initials for manage spot size in the plots
var_height <- 24 # 24/3=8
var_width <- 26 # 36/4=9
var_point_size <- 1.5

set.seed(07112024)

lst_order <- sort(unique(spe$sample_id))

vis_grid_clus(
  spe = spe,
  clustervar = "BayesSpace",
  sample_order = lst_order,
  height = var_height, # 8
  width = var_width, # 9
  point_size = var_point_size,
  pdf = here(dir_plots, "spe_BayesSpace_annotated_grid.pdf"),
  sort_clust = FALSE,
  guide_point_size = 10,
  colors = c("#b2df8a", "#e41a1c", "#377eb8", "#4daf4a", "#ff7f00", "gold", "#a65628",
             "#999999", "black", "grey", "purple", "white", "red",
  )
  
  # vis_clus(
  #   spe,
  #   sampleid = unique(spe$sample_id)[1],
  #   clustervar = "BayesSpace",
  #   colors = c("#b2df8a", "#e41a1c", "#377eb8", "#4daf4a", "#ff7f00", "gold", "#a65628",
  #              "#999999", "black", "grey", "purple", "white", "red"),
  #   spatial = TRUE,
  #   image_id = "lowres",
  #   alpha = NA,
  #   point_size = var_point_size,
  #   pdf = here(dir_plots, "spe_BayesSpace_annotated_individual_sample.pdf"),
  #   auto_crop = TRUE,
  #   na_color = "#CCCCCC40",
  #   is_stitched = FALSE,
  #   guide_point_size = 10
  # )
  
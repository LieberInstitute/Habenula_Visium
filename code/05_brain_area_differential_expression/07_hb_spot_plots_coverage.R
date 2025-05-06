library("here")
library("spatialLIBD") #[1] ‘1.21.4’
library("purrr")
library("sessioninfo")

## load spe data

spe_in <- here("processed-data", "04_harmony_BayesSpace", "spe_harmony_ann.rds")
spe <- readRDS(spe_in)


## Quick inspection
spe
# colnames(colData(spe))
levels(colData(spe)$BayesSpace)
# [1] "Sp13D01.Oligo"      "Sp13D02.Endo"       "Sp13D03.Endo"      
# [4] "Sp13D04.Astrocyte"  "Sp13D05.Inhib.Thal" "Sp13D06.Astrocyte" 
# [7] "Sp13D07.Astrocyte"  "Sp13D08.OPC_Astroc" "Sp13D09.Astrocyte_"
# [10] "Sp13D10.Oligo"      "Sp13D11.Habenula"   "Sp13D12.Oligo"     
# [13] "Sp13D13.Oligo"  
## Get number of TRUE spots in tissue
map(unique(spe$sample_id), ~ summary(spe$in_tissue[spe$sample_id == .x] == TRUE))
unique(spe$sample_id)

## set dir input 

dir_plots <- here("plots","05_brain_area_differential_expression", "07_hb_spot_plots_coverage")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

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
  guide_point_size = 5,
  colors = c("grey", "#b2df8a", "#e41a1c", "#377eb8", "#4daf4a",
             "#ff7f00", "black", "#a65628", "#999999", "blue",
             "purple", "gold")
  )
  
vis_grid_clus(
  spe = spe,
  clustervar = "BayesSpace",
  sample_order = lst_order,
  height = var_height, # 8
  width = var_width, # 9
  point_size = var_point_size,
  pdf = here(dir_plots, "spe_BayesSpace_annotated_hb_grid.pdf"),
  sort_clust = FALSE,
  guide_point_size = 5,
  colors = c("grey", "grey", "grey", "grey", "grey",
             "grey", "grey", "grey", "grey", "blue",
             "grey", "grey")
)





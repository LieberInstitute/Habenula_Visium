library("spatialLIBD")
library("here")
library("tidyverse")
# library("scran")
# library("BiocParallel")
# library("scater")
# library("scry")
# library("BiocSingular")
library("sessioninfo")
#library("HDF5Array")


dir_rdata <- here("processed-data", "04_harmony_BayesSpace")
filtered_in_path <- file.path(dir_rdata, "spe_qc_filtered_logcounts.rds")
#filtered_hdf5_dir <- file.path(dir_rdata, "spe_filtered_hdf5")
dir_plots <- here("plots", "04_harmony_BayesSpace", "one_gene_specific_kristen")

## Create output directories
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

## load a filtered spe object
spe <- readRDS(filtered_in_path)

cat("Number of spots after removed any remaining empty spots and/or genes with zero counts:", dim(spe)[2],"\n")

colnames(colData(spe))

# lst_viridis <- c(gene_name_ID_OPRM1 = "Gene_OPRM1_logCounts-viridis.pdf", 
#                  gene_name_ID_POU4F1 = "Gene_POU4F1_logCounts-viridis.pdf")
# 
# lst_magma <- c(OPRM1 = "Gene_OPRM1_logCounts-magma.pdf", 
#                  POU4F1 = "Gene_POU4F1_logCounts-magma.pdf")
# 
# lst_viridis_FALSE <- c(OPRM1 = "Gene_OPRM1_logCounts-spatialFALSE-viridis.pdf",
#                        POU4F1 = "Gene_POU4F1_logCounts-spatialFALSE-viridis.pdf")
# 
# lst_magma_FALSE <- c(OPRM1 = "Gene_OPRM1_logCounts-spatialFALSE-magma.pdf", 
#                      POU4F1 = "Gene_POU4F1_logCounts-spatialFALSE-magma.pdf")


# map2(as.vector(names(lst_viridis)), as.vector(lst_viridis), ~ vis_grid_gene(
#   spe = spe,
#   geneid = .x,
#   point_size = 4,
#   cont_colors = viridisLite::magma(21, direction = -1),
#   spatial = TRUE,
#   pdf = here(dir_plots, .y)
# ))

lst_Habenula <- rowData(spe)$gene_search[
  rowData(spe)$gene_name %in% 'OPRM1'
]

pdf_name <- "Gene_OPRM1_logCounts-viridis_reverse.pdf" 

vis_grid_gene(
  spe = spe,
  geneid = lst_Habenula,
  point_size = 4,
  cont_colors = viridisLite::viridis(21, direction = 1),
  spatial = TRUE,
  pdf = here(dir_plots, pdf_name)
)

lst_Habenula <- rowData(spe)$gene_search[
  rowData(spe)$gene_name %in% 'POU4F1'
]

pdf_name <- "Gene_POU4F1_logCounts-viridis_reverse.pdf" 

vis_grid_gene(
  spe = spe,
  geneid = lst_Habenula,
  point_size = 4,
  cont_colors = viridisLite::viridis(21, direction = 1),
  spatial = TRUE,
  pdf = here(dir_plots, pdf_name),
)



# vis_grid_gene(
#   spe = spe,
#   geneid = lst_Habenula,
#   point_size = 4,
#   cont_colors = viridisLite::turbo(21, direction = -1),
#   spatial = TRUE,
#   pdf = here(dir_plots, pdf_name),
#   assayname = counts
# )



######################################################################################################




## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

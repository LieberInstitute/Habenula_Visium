library(here)
library(SpatialExperiment)
library(Seurat)
library(HDF5Array)
library(nnSVG)
library(PRECAST)
library(tidyverse)
library(Matrix)
library(sessioninfo)


## Run PRECAST for integrating and analyzing multiple spatially resolved transcriptomics (SRT) datasets. It unifies spatial factor analysis simultaneously with spatial clustering and embedding alignment, requiring only partially shared cell/domain clusters across datasets.

k = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# test k <- 3

# PRECAST only requests the spe counts assay (with tissue spot QCed / spe_qcED_spatialLIBD.rds). It can be RDS or hdf5 object.
spe_dir <- here("processed-data", "04_harmony_BayesSpace", "spe_qcED_spatialLIBD_log.rds")  # "spe_qcED_spatialLIBD.rds" 
out_path = here('processed-data', '08_precast')
dir.create(spe_dir, showWarnings = FALSE)
dir.create(out_path, showWarnings = FALSE)
# num_genes = 2000

spe = readRDS(spe_dir)  # loadHDF5SummarizedExperiment(spe_dir)
# spe
# class: SpatialExperiment 
# dim: 25826 16613 
# metadata(0):
#   assays(2): counts logcounts
# names(colData(spe))

spe$row = spe$array_row
spe$col = spe$array_col
## PRECAST expects array coordinates in 'row' and 'col' columns. You should use transformed col/row coordinates  if you have stitched data
# spe$row = spe$array_row_transformed
# spe$col = spe$array_col_transformed
## This is other alternative with no stitched data  
## https://github.com/LieberInstitute/spatialDLPFC/blob/bd93c980d7653579f81ff1c91c309cea0c7474a6/code/analysis/01_build_spe/01_build_spe.R#L653-L657

## Create a list of Seurat objects: one per sample_id 
seu_list = lapply(
  unique(spe$sample_id),
  function(sample) {  
    small_spe = spe[, spe$sample_id == sample]
  
    CreateSeuratObject(
      #   Bring into memory to greatly improve speed
      counts = as(assays(small_spe)$counts, "dgCMatrix"),
      meta.data = as.data.frame(colData(small_spe)),
      project = 'Hb_Visium'
    )
  }
)


################################################################################
## This chunk is for is for running precast with HVDG (High Variable Deviance Genes) 
##      obtained from GLM-PCA (Binomial model .20p (~2000 genes - count assay)
################################################################################

## Input our previous selected HVDG (High Variable Deviance Genes) processed to 
##      compute GLM-PCA (Binomial model 2000 genes - count assay)
dir_rdata <- here("processed-data", "04_harmony_BayesSpace", "hdgs.hb.Rdata")
hvdg_all <- load(dir_rdata)
hvdg <- hdgs.hb.2000


################################################################################
## This chunk is for running precast with SVG (Spatial Variable Genes) 
##      SVG are calculated with nnSVG on pre-calculated log assays with scran
## MOVE THIS CHUNK TO AN INDEPENDENT SCRIPT TO AVOID RE-CALCULATE THE SVGs  
################################################################################

# ## Run nnSVG() function as we have log normalized count in our object.
# #spe <- nnSVG(spe, X = model.matrix(~capture_area_))
# spe <- nnSVG(spe)
# 
# # Quick inspection to results
# rowData(spe)
# # number of significant SVGs
# table(rowData(spe)$padj <= 0.05)
# # show results for top n SVGs
# n <- 100
# svgs <- rowData(spe)[order(rowData(spe)$rank)[1:n], ]
# print(paste("Top 100 SVGs at pajd < 0.05:", length(svgs)))
# 
# saveRDS(svgsfile = file.path(dir_rdata, "top.svg.RDS")





## Prepare the PRECASTObject with preprocessing step based on the Seurat list object seuList.
##  1. Filter low-quality spots and genes
##  2. Select the top 2000 variable genes (by setting gene.number=2000) for each data batch using FindVariableFeatures() function in Seurat package for highly variable genes.
##  3. Conduct strict quality control for data_filter2 by filtering spots and genes, controlled by the arguments .

# Note. If the argument customGenelist is not NULL, then this function only does (3) based on customGenelist gene list.
#row.names(seu_list[[1]])

set.seed(31072024)

# Note. Set the correct assay as Default. 

pre_obj = CreatePRECASTObject(
  seuList = seu_list,
  selectGenesMethod = NULL,
  customGenelist = hvdg
  #customGenelist = svgs
)

pre_obj@seulist

## Add adjacency matrix list for a PRECASTObj object to prepare for PRECAST model fitting.

#   Setting platform to "Visium" just means to use array indices, which should
#   work fine despite the abnormal/ "artificial" capture area we've created by
#   stitching
pre_obj <- AddAdjList(pre_obj, platform = "Visium")

#   Following https://feiyoung.github.io/PRECAST/articles/PRECAST.BreastCancer.html,
#   which involves overriding some default values, though the implications are not
#   documented
pre_obj <- AddParSetting(
  pre_obj, Sigma_equal = FALSE, verbose = TRUE, maxIter = 30
)
pre_obj@parameterList

#   Fit model. users can specify the number of clusters 𝐾 or input an integer vector. e.g. K=6:9
pre_obj <- PRECAST(pre_obj, K = k)
# -----Intergrative data info.: 5 samples, 1998 genes X 16604 spots------
#   -----PRECAST model setting: error_heter=TRUE, Sigma_equal=FALSE, Sigma_diag=TRUE, mix_prop_heter=TRUE
# Start computing intial values... 

## Select a best model and re-organize the results
pre_obj <- SelectModel(pre_obj)
## backup the fitted results in resList
resList <- PRECASTObj@resList
str(PRECASTObj@resList)
## Integrate multiple SRT data based on the PreCast object by PreCast modeling fitting
##    At this point we should probably look at other biological covariates to oncisder when removing batch effect: 
# PRECASTObj@seulist[[1]]@meta.data / covariates_use=NULL
pre_obj = IntegrateSpaData(pre_obj, species = "Human")


##   Extract PRECAST results, clean up column names, and export to CSV
precast_results <- pre_obj@meta.data |>
    rownames_to_column("key") |>
    as_tibble() |>
    select(-orig.ident) |>
    rename_with(~ sub('_PRE_CAST', '', .x)) |>
write_csv(precast_results, file.path(out_path, sprintf("PRECAST_k%s.csv", k))) 

  
saveRDS(pre_obj, file.path(out_path, sprintf("/PRECAST_k%s_integrated.rds", k)))
out_path <- dirname(out_path)
saveRDS(resList, file.path(out_path, sprintf("/restlist_k%s_integrated.rds", k)))

session_info()



################################################################################
## MOVE THIS CHUNK TO AN INDEPENDENT SCRIPT TO ONLY BUID PLOTS FOR EXPLORATION
################################################################################

plot_dir = here('plots', '08_precast')

## Some visualizations
cols_cluster <- chooseColors(palettes_name = 'Nature 10', n_colors = 7, plot_colors = TRUE)
plt_precast <- SpaPlot(pre_obj, batch=NULL, cols=cols_cluster, point_size=2, combine=TRUE)
plt_precast

pre_obj <- RunTSNE(pre_obj, reduction = "PRECAST", tSNE.method = "FIt-SNE")

p1 <- dimPlot(pre_obj, item = "cluster", point_size = 0.5, font_family = "serif", cols = cols_cluster,
              border_col = "gray10", nrow.legend = 14, legend_pos = "right")
# 
# cols_batch <- chooseColors(palettes_name = "Classic 20", n_colors = 10, plot_colors = TRUE)
# p2 <- dimPlot(pre_obj, item = "batch", point_size = 0.5, font_family = "serif", cols = cols_batch,
#               border_col = "gray10", nrow.legend = 14, legend_pos = "right")
# 
# pdf(file.path(plot_dir, "/precast_tSNE.pdf"), width = 12, height = 8)
# plot_grid(p1, p2, ncol = 2)
# dev.off()
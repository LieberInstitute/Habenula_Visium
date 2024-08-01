library(here)
library(PRECAST)
library(HDF5Array)
library(Seurat)
library(sessioninfo)
library(tidyverse)
library(Matrix)
library(SpatialExperiment)

## Run PRECAST for integrating and analyzing multiple spatially resolved transcriptomics (SRT) datasets. It unifies spatial factor analysis simultaneously with spatial clustering and embedding alignment, requiring only partially shared cell/domain clusters across datasets.

k = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# test k <- 3

# spe_dir = here(
#   'processed-data', '05_harmony_BayesSpace', 'spe_filtered_hdf5'
# )
# spe with in-tissue spot counts, QCed and excluding those manually annotated as low quality or having tissue artifacts.
spe_dir <- here("processed-data", "04_harmony_BayesSpace", "spe_qcED_spatialLIBD.rds")

out_path = here('processed-data', '08_precast', paste0('PRECAST_k', k, '.csv'))
# num_genes = 2000

set.seed(31072024)
dir.create(dirname(out_path), showWarnings = FALSE)

#spe = loadHDF5SummarizedExperiment(spe_dir)
spe = readRDS(spe_dir)
# spe
# class: SpatialExperiment 
# dim: 25826 16613 
# metadata(0):
#   assays(1): counts

names(colData(spe))

#   PRECAST expects array coordinates in 'row' and 'col' columns
# spe$row = spe$array_row_transformed
# spe$col = spe$array_col_transformed

# note csc. need double chk if I need to transform the coordinates 
spe$row = spe$array_row
spe$col = spe$array_col

## Alternative 2
## https://github.com/LieberInstitute/spatialDLPFC/blob/bd93c980d7653579f81ff1c91c309cea0c7474a6/code/analysis/01_build_spe/01_build_spe.R#L653-L657
# auto_offset_row <- as.numeric(factor(unique(spe$sample_id))) * 100
# names(auto_offset_row) <- unique(spe$sample_id)
# spe$row <- colData(spe)$array_row + auto_offset_row[spe$sample_id]
# spe$col <- colData(spe)$array_col

#   Create a list of Seurat objects: one per sample_id 
seu_list = lapply(
  # unique(spe$donor),
  unique(spe$sample_id),
  # function(donor) {
    # small_spe = spe[, spe$donor == donor]
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


## input our previous selected HVDG (High Variable Deviance Genes) processed to 
##      compute GLM-PCA (Binomial model 2000 genes - count assay)
dir_rdata <- here("processed-data", "04_harmony_BayesSpace", "hdgs.hb.Rdata")
hvdg_all <- load(dir_rdata)
hvdg <- hdgs.hb.2000

# svgs = read.csv(svg_path) |>
#   as_tibble() |>
#   arrange(nnsvg_avg_rank_rank) |>
#   slice_head(n = num_genes) |>
#   pull(gene_id)

## Prepare the PRECASTObject with preprocessing step based on the Seurat list object seuList.
##  1. Filter low-quality spots and genes
##  2. Select the top 2000 variable genes (by setting gene.number=2000) for each data batch using FindVariableFeatures() function in Seurat package for highly variable genes.
##  3. Conduct strict quality control for data_filter2 by filtering spots and genes, controlled by the arguments .

# Note. If the argument customGenelist is not NULL, then this function only does (3) based on customGenelist gene list.
#row.names(seu_list[[1]])

pre_obj = CreatePRECASTObject(
  seuList = seu_list,
  selectGenesMethod = NULL,
  customGenelist = hvdg
  #customGenelist = svgs
  

  #   Using defaults for gene-filtering-related parameters. Though each donor
  #   consists of more spots than 1 typical Visium capture area (and would
  #   thus be expected to throw off the appropriateness of the defaults for
  #   'premin.spots', etc), we're using SVGs from nnSVG as input, and these
  #   genes already passed a similar reasonable expression cutoff:
  #   https://github.com/LieberInstitute/spatial_NAc/blob/61d1e198536a80bddca93017ea6eb8169af5d978/code/05_harmony_BayesSpace/05-run_nnSVG.R#L40-L45
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

#   Fit model. users can specify the number of clusters 𝐾
pre_obj <- PRECAST(pre_obj, K = k)
# -----Intergrative data info.: 5 samples, 1998 genes X 16604 spots------
#   -----PRECAST model setting: error_heter=TRUE, Sigma_equal=FALSE, Sigma_diag=TRUE, mix_prop_heter=TRUE
# Start computing intial values... 

## Select a best model and re-organize the results
pre_obj <- SelectModel(pre_obj)
str(PRECASTObj@resList)
## Integrate data
pre_obj = IntegrateSpaData(pre_obj, species = "Human")

## Some visualizations
cols_cluster <- chooseColors(palettes_name = 'Nature 10', n_colors = 7, plot_colors = TRUE)
plt_precast <- SpaPlot(pre_obj, batch=NULL, cols=cols_cluster, point_size=2, combine=TRUE)
plt_precast

#   Extract PRECAST results, clean up column names, and export to CSV
pre_obj@meta.data |>
  rownames_to_column("key") |>
  as_tibble() |>
  select(-orig.ident) |>
  rename_with(~ sub('_PRE_CAST', '', .x)) |>
  write_csv(out_path)

session_info()

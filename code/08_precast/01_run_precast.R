library(here)
library(PRECAST)
library(HDF5Array)
library(Seurat)
library(sessioninfo)
library(tidyverse)
library(Matrix)
library(SpatialExperiment)

k = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

spe_dir = here(
  'processed-data', '05_harmony_BayesSpace', 'spe_filtered_hdf5'
)
svg_path = here(
  'processed-data', '05_harmony_BayesSpace', 'nnSVG_out',
  'summary_across_samples.csv'
)
out_path = here('processed-data', '10_precast', paste0('PRECAST_k', k, '.csv'))
num_genes = 2000

set.seed(1)
dir.create(dirname(out_path), showWarnings = FALSE)

spe = loadHDF5SummarizedExperiment(spe_dir)

#   PRECAST expects array coordinates in 'row' and 'col' columns
spe$row = spe$array_row_transformed
spe$col = spe$array_col_transformed

#   Create a list of Seurat objects: one per donor
seu_list = lapply(
  unique(spe$donor),
  function(donor) {
    small_spe = spe[, spe$donor == donor]
    
    CreateSeuratObject(
      #   Bring into memory to greatly improve speed
      counts = as(assays(small_spe)$counts, "dgCMatrix"),
      meta.data = as.data.frame(colData(small_spe)),
      project = 'spatialNAc'
    )
  }
)

svgs = read.csv(svg_path) |>
  as_tibble() |>
  arrange(nnsvg_avg_rank_rank) |>
  slice_head(n = num_genes) |>
  pull(gene_id)

pre_obj = CreatePRECASTObject(
  seuList = seu_list,
  selectGenesMethod = NULL,
  customGenelist = svgs
  #   Using defaults for gene-filtering-related parameters. Though each donor
  #   consists of more spots than 1 typical Visium capture area (and would
  #   thus be expected to throw off the appropriateness of the defaults for
  #   'premin.spots', etc), we're using SVGs from nnSVG as input, and these
  #   genes already passed a similar reasonable expression cutoff:
  #   https://github.com/LieberInstitute/spatial_NAc/blob/61d1e198536a80bddca93017ea6eb8169af5d978/code/05_harmony_BayesSpace/05-run_nnSVG.R#L40-L45
)

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

#   Fit model
pre_obj <- PRECAST(pre_obj, K = k)
pre_obj <- SelectModel(pre_obj)
pre_obj = IntegrateSpaData(pre_obj, species = "Human")

#   Extract PRECAST results, clean up column names, and export to CSV
pre_obj@meta.data |>
  rownames_to_column("key") |>
  as_tibble() |>
  select(-orig.ident) |>
  rename_with(~ sub('_PRE_CAST', '', .x)) |>
  write_csv(out_path)

session_info()

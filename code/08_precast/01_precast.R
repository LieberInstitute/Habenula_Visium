library(getopt)
library(sessioninfo)
library(here)
library(PRECAST)
library(HDF5Array)
library(Seurat)
library(tidyverse)
library(Matrix)
library(SpatialExperiment)

# Import command-line parameters
spec <- matrix(
  c(
    c("final_step", "k"),
    c("f", "k"),
    rep("1", 2),
    c("character", "numeric"),
    c(
      "Using array coordinates after either 'imagej' or 'samui",
      "Number of clusters"
    )
  ),
  ncol = 5
)
opt <- getopt(spec)

print("Using the following parameters:")
print(opt)

spe_dir = here('processed-data', '15_samui_imagej_comparison', 'spe')
svg_path = here(
  'processed-data', '05_harmony_BayesSpace', '07-run_nnSVG', 'nnSVG_out',
  'summary_across_samples.csv'
)
out_path = here(
  'processed-data', '15_samui_imagej_comparison', 'precast_out',
  sprintf('PRECAST_k%s_%s.csv', opt$k, opt$final_step)
)
num_genes = 2000

set.seed(1)
dir.create(dirname(out_path), showWarnings = FALSE)

spe = loadHDF5SummarizedExperiment(spe_dir)

#   PRECAST expects array coordinates in 'row' and 'col' columns
spe$row = spe[[paste0('array_row_', opt$final_step)]]
spe$col = spe[[paste0('array_col_', opt$final_step)]]

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
pre_obj <- PRECAST(pre_obj, K = opt$k)
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

## This script was made using slurmjobs version 1.2.2
## available from http://research.libd.org/slurmjobs/
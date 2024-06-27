# library(slurmjobs)
# slurmjobs::job_single('02_explore_expr_variability', create_shell = TRUE, memory = '20G', command = "02_explore_expr_variability.R")

# To submit the job use: sbatch 02_explore_expr_variability.sh

k <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## For testing
if (FALSE) {
  k <- 2
}

library("here")
library("sessioninfo")
library("SingleCellExperiment")
library("scater")
## Load BayesSpace colors
library(Polychrome)

colors_bayesSpace <- Polychrome::palette36.colors(28)
names(colors_bayesSpace) <- c(1:28)

#source(here("code", "analysis", "colors_bayesSpace.R"), echo = TRUE, max.deparse.length = 500)
names(colors_bayesSpace) <-
  paste0("Sp", sprintf("%02d", k), "D", sprintf("%02d", as.integer(names(colors_bayesSpace))))

## output directory
dir_rdata <- here("processed-data", "05_layer_differential_expression")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dir_rdata)) ## Check that it was created successfully
dir_plots <- here("plots", "05_layer_differential_expression")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dir_plots))

## load sce_pseudo data
sce_pseudo <-
  readRDS(
    file.path(
      dir_rdata,
      paste0("sce_pseudo_BayesSpace_k", sprintf("%02d", k), ".rds")
    )
  )

## Define variables to use
vars <- c(
  "age",
  "sample_id",
  "BayesSpace",
  "subject",
  "sex" 
)

## Plot PCs with different colors
## Each point here is a sample
reducedDim(sce_pseudo)

pdf(file = file.path(dir_plots, paste0("sce_pseudo_PCs_k", sprintf("%02d", k), ".pdf")), width = 14, height = 14)
for (var in vars) {
  p <- plotPCA(
    sce_pseudo,
    colour_by = var,
    ncomponents = 12,
    point_size = 1,
    label_format = c("%s %02i", " (%i%%)"),
    percentVar = metadata(sce_pseudo)$PCA_var_explained
  )
  if (var == "BayesSpace") {
    p <- p + scale_color_manual("BayesSpace", values = colors_bayesSpace)
  }
  print(p)
}
dev.off()


## Obtain percent of variance explained at the gene level
## using scater::getVarianceExplained()
variance_expl <- getVarianceExplained(sce_pseudo,
                             variables = vars
)
head(variance_expl)
summary(variance_expl)

## Now visualize the percent of variance explained across all genes
pdf(file = file.path(dir_plots, paste0("sce_pseudo_gene_explanatory_vars_k", sprintf("%02d", k), ".pdf")))
plotExplanatoryVariables(variance_expl)
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

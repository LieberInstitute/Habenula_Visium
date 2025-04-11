library("spatialLIBD")
library("markdown")
library("here")

here::here("code", "03_spatialLIBD_app_pseudobulk")

## spatialLIBD uses golem.
## Golem is a framework for building production-grade shiny applications
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

######### Load All required objects #########

## Set BayesSpace k selection
BayesSpace_k <- "09"

# ## load harmony_BayesSpace spe object
spe <- readRDS("spe_pseudobulk_shiny.rds")
# lobstr::obj_size(spe)
# 4.03 GB

colnames(colData(spe))
#head(spe$BayesSpace_harmony_k16)
# NULL

## load the pseudobulk object sce_pseudo
# sce_pseudo_name <- paste0("sce_pseudo_BayesSpace_k", BayesSpace_k,".rds")
sce_pseudo_name <- paste0("sce_pseudo_PCA_brain_area_k", BayesSpace_k,".rds")
sce_pseudo <- readRDS(sce_pseudo_name)

## load modeling results for any k9 clustering/pseudobulking
modeling_results_name <- paste0("modeling_results_BayesSpace_k", BayesSpace_k ,".Rdata")
load(modeling_results_name, verbose = TRUE)

## Load all significant genes
signif_genes_name <- paste0("sig_genes_k", BayesSpace_k,".Rdata")
load(signif_genes_name, verbose = TRUE)

# Quick inspection
colnames(colData(spe))

## Import `BayesSpace_harmony_k%` clusters and assign them to `BayesSpace` new column
spe <- cluster_import(spe,
    cluster_dir = "clusters_BayesSpace",
    prefix = ""
)
# Overwriting 'spe$key'. Set 'overwrite = FALSE' if you do not want to overwrite it

BS_k_column <- paste0("BayesSpace_harmony_k", BayesSpace_k)
spe$BayesSpace <- spe[[BS_k_column]]
#colnames(colData(spe))
#head(spe$BayesSpace)

## Quickly explore the data
vars <- colnames(colData(spe))
colnames(colData(spe)) <- vars <- gsub("X10x", "10x", vars)

colors_BayesSpace <- Polychrome::palette36.colors(28)
names(colors_BayesSpace) <- c(1:28)
m <- match(as.character(spe$BayesSpace_harmony_k09), names(colors_BayesSpace))
stopifnot(all(!is.na(m)))
spe$BayesSpace_colors <- spe$BayesSpace_harmony_k09_colors <- colors_BayesSpace[m]


title_name <- paste0("spatialHabenula, Visium, Sp", BayesSpace_k)

spatialLIBD::run_app(
    spe,
    sce_layer = sce_pseudo,
    modeling_results = modeling_results,
    sig_genes = sig_genes,
    title = title_name,
    spe_discrete_vars = c( ## this are the variables for the spe object not the sce_pseudo object
        "BayesSpace",
        "ManualAnnotation",
        "overlaps_tissue",
        vars[grep("^10x_", vars)],
        vars[grep("^scran_", vars)],
        "edge_spot",
        vars[grep("^SNN_k10", vars)],
        vars[grep("^BayesSpace_harmony_", vars)],
        "BayesSpace_colors"
    ),
    spe_continuous_vars = c(
        "sum_umi",
        "sum_gene",
        "expr_chrM",
        "expr_chrM_ratio",
        "edge_distance"
    ),
    default_cluster = "BayesSpace",
    docs_path = "www"
)


## Note. If fails to read the rds object, go to Session Menu -> Set Working Directory -> To source File location

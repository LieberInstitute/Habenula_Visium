library(spatialLIBD)
library(markdown)
library(here)
library(tidyverse)
library(qs2)

#   For interactive testing at JHPCE
# setwd(here('code', '12_apps_and_sharing', 'shiny_app'))

spe_path = 'spe_shiny.qs2'
modeling_path = 'modeling_results.rds'
spe_pb_path = 'spe_pb_shiny.qs2'
sig_genes_path = 'sig_genes_shiny.qs2'
docs_dir = 'www'
discrete_vars = c(
    'sample_id', 'donor', 'tissue_piece', 'ManualAnnotation', 'ficture_cluster',
    'banksy_cluster', 'cell_type'
)
continuous_vars = c(
    'bin_count', 'sum_umi', 'sum_gene', 'expr_chrM', 'expr_chrM_ratio'
)

## spatialLIBD uses golem
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

#   Load objects
spe = qs_read(spe_path)
modeling_results = readRDS(modeling_path)
spe_pb = readRDS(spe_pb_path)
sig_genes = qs_read(sig_genes_path)

## Deploy the website
run_app(
    spe,
    sce_layer = spe_pb,
    modeling_results = modeling_results,
    sig_genes = sig_genes,
    title = "habenula_atlas_Visium_HD",
    spe_discrete_vars = discrete_vars,
    spe_continuous_vars = continuous_vars,
    default_cluster = "cell_type",
    docs_path = docs_dir,
    is_stitched = TRUE
)

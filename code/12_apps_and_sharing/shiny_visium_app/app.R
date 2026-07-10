library(spatialLIBD)
library(markdown)
library(here)
library(tidyverse)
library(qs2)

#   For interactive testing at JHPCE
# setwd(here('code', '12_apps_and_sharing', 'shiny_visium_app'))

discrete_vars = c(
    "ManualAnnotation", "brain_id", "sex", "ethnicity",
    sprintf("SNN_k10_k%d", seq(4, 28)),
    sprintf("BayesSpace_harmony_k%02d", seq(2, 28))
)
continuous_vars = c(
    "sum_umi", "sum_gene", "expr_chrM", "expr_chrM_ratio", "age", "pmi", "rin"
)

## spatialLIBD uses golem
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

#   Load objects
spe = qs_read('spe_shiny.qs2')
modeling_results = get(load('modeling_results.Rdata'))
sce_pb = qs_read('sce_pb_shiny.qs2')
sig_genes = qs_read('sig_genes_shiny.qs2')

## Deploy the website
run_app(
    spe,
    sce_layer = sce_pb,
    modeling_results = modeling_results,
    sig_genes = sig_genes,
    title = "habenula_atlas_Visium",
    spe_discrete_vars = discrete_vars,
    spe_continuous_vars = continuous_vars,
    default_cluster = "BayesSpace_harmony_k09",
    docs_path = 'www'
)

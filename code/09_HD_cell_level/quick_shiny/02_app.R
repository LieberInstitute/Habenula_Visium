library(spatialLIBD)
library(markdown)
library(here)
library(tidyverse)
library(qs2)

spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'quick_shiny',
    'spe_shiny.qs2'
)
docs_dir = here('code', '09_HD_cell_level', 'quick_shiny', 'www')
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

spe = qs_read(spe_path)

## Deploy the website
run_app(
    spe,
    sce_layer = NULL,
    modeling_results = NULL,
    sig_genes = NULL,
    title = "habenula_atlas_Visium_HD",
    spe_discrete_vars = discrete_vars,
    spe_continuous_vars = continuous_vars,
    default_cluster = "cell_type",
    docs_path = docs_dir,
    is_stitched = TRUE
)

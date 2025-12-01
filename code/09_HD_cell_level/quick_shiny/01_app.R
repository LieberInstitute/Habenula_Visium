library(spatialLIBD)
library(markdown)
library(here)
library(tidyverse)

spe_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'spe_norm_filtered.rds'
)
docs_dir = here('code', '09_HD_cell_level', 'quick_shiny', 'www')
banksy_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'banksy', 'lambda0_2',
    'leiden_res0_1.csv'
)

## spatialLIBD uses golem
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

spe <- readRDS(spe_path)

spe$banksy = tibble(key = spe$key) |>
    left_join(read_csv(banksy_path, show_col_types = FALSE), by = "key") |>
    pull(banksy_lambda0_2)
stopifnot(!any(is.na(spe$banksy)))

## Deploy the website
run_app(
    spe,
    sce_layer = NULL,
    modeling_results = NULL,
    sig_genes = NULL,
    title = "habenula_atlas_HD",
    spe_discrete_vars = c(
        "ManualAnnotation",
        "labels_joint_source",
        "banksy"
    ),
    spe_continuous_vars = c(
        "sum_umi",
        "sum_gene",
        "expr_chrM",
        "expr_chrM_ratio"
    ),
    default_cluster = "banksy",
    docs_path = docs_dir,
    is_stitched = TRUE
)

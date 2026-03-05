library(spatialLIBD)
library(markdown)
library(tidyverse)
library(here)

#   At JHPCE
setwd(here("code", "12_HD_shiny"))

options("golem.app.prod" = TRUE)
options(repos = BiocManager::repositories())

spe = readRDS("spe_norm_filtered_split.rds")

anno_df = read_csv("cluster_annotation.csv", show_col_types = FALSE)

#   Tibble including banksy cluster and cell type
extra_coldata = tibble(key = spe$key) |>
    left_join(
        read_csv("leiden_res1_8.csv", show_col_types = FALSE), by = 'key'
    ) |>
    mutate(
        cell_type = anno_df$fine_cell_type[
            match(as.character(banksy), anno_df$cluster)
        ]
    )
stopifnot(!any(is.na(extra_coldata$cell_type)))

spe$banksy = extra_coldata$banksy
spe$cell_type = extra_coldata$cell_type

# lobstr::obj_size(spe) is 3.95GB as of 2026-03-05

vars <- colnames(colData(spe))
run_app(
    spe = spe,
    sce_layer = NULL,
    modeling_results = NULL,
    sig_genes = NULL,
    spe_discrete_vars = c(
        "ManualAnnotation",
        "exclude_overlapping",
        "banksy",
        "cell_type"
    ),
    spe_continuous_vars = c(
        "bin_count",
        "sum_umi",
        "sum_gene",
        "expr_chrM",
        "expr_chrM_ratio"
    ),
    default_cluster = "cell_type",
    docs_path = "www",
    is_stitched = TRUE
)

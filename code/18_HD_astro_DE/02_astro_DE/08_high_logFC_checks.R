library(tidyverse)
library(here)
library(sessioninfo)
library(qs2)
library(spatialLIBD)
library(edgeR)
library(limma)
library(duckplyr)

de_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'astro_DE',
    'main_results', 'DE_0.parquet'
)
dge_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'astro_DE',
    'main_results', 'DGE.qs2'
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dge = qs_read(dge_path)
de_df = read_parquet_duckdb(de_path) |>
    collect()

de_highest_df = de_df |>
    arrange(desc(logFC)) |>
    slice_head(n = 1)
gene_id = de_highest_df$gene_id

y <- dge$EList$E[gene_id, ]
w <- dge$EList$weights[match(gene_id, rownames(dge$EList$E)), ]

manual_fit <- lm(
    y ~ astro_label + ncells + expr_chrM_ratio,
    data = dge$targets,
    weights = w
)

all.equal(
    dge$coefficients[gene_id, 'astro_labellateral'],
    unname(coef(manual_fit)['astro_labellateral']),
    de_highest_df$logFC
)

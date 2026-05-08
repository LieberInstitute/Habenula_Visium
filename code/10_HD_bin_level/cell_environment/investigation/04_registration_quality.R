#   It looks like the issue is related to quality of the extracellular
#   registration stats. There are nearly no significant genes for each
#   cluster, and any correlation heatmaps using the stats have tiny
#   correlations spatially unrelated to where the clusters belong. Do
#   FICTURE's reported top genes relate at all to top genes from the
#   registration process?
#
#   I'm fairly sure there won't be a substantial relationship. Quickly run
#   some tests with AI

library(here)
library(sessioninfo)
library(tidyverse)

#   Interactively tried every k actually
k = 10

ficture_genes_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'ficture_outputs', 'cleaningy', sprintf('k_%d', k), 'analysis',
    sprintf('nF%d.d_12', k),
    sprintf('nF%d.d_12.decode.prj_12.r_4_5.factor.info.tsv', k)
)
model_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'registration', 'modeling_results', sprintf('%s.rds', k)
)

ficture_genes = read_tsv(ficture_genes_path)
model = readRDS(model_path)
enrich = model$enrichment

#   For each FICTURE factor, test whether the registration t-stats agree with
#   FICTURE's gene weight ranking using Spearman rank correlation.
#   If both methods capture real signal, genes ranked highly by FICTURE weight
#   should also have high t-stats in the matched registration cluster (positive rho).
#   Near-zero or negative rho indicates the registration contains no convergent signal.

ficture_top = ficture_genes |>
    pivot_longer(
        cols = starts_with('TopGene'),
        names_to = 'ranking_method',
        names_pattern = 'TopGene_(.*)',
        values_to = 'genes'
    )

spearman_results = ficture_top |>
    filter(ranking_method == 'weight', genes %in% rownames(enrich)) |>
    mutate(ficture_rank = row_number(), .by = Factor) |>
    mutate(t_col = paste0('t_stat_X', Factor)) |>
    rowwise() |>
    mutate(t_stat = enrich[genes, t_col]) |>
    ungroup() |>
    summarise(
        n = n(),
        spearman_r = cor(ficture_rank, t_stat, method = 'spearman'),
        .by = Factor
    ) |>
    arrange(Factor)

#   10X announced that many probes in the probe set as input to SpaceRanger were
#   "bad". This script investigates the prevalence of bad probes in top markers
#   for Banksy and FICTURE clusters

library(here)
library(tidyverse)
library(spatialLIBD)
library(HDF5Array)
library(sessioninfo)

spe_dir = here('processed-data', '09_HD_cell_level', 'spe_norm_filtered')
banksy_model_path = here(
    'processed-data', '09_HD_cell_level', 'registration_banksy',
    'modeling_results', '1.rds'
)
ficture_model_path = here(
    'processed-data', '10_HD_bin_level', 'ficture_harmony',
    'modeling_results', 'library_normalized.rds'
)
bad_path = here(
    'processed-data', '10_HD_bin_level', 'bad_probe_genes',
    'OtherGenesWithExcludedProbes.txt'
)
very_bad_path = here(
    'processed-data', '10_HD_bin_level', 'bad_probe_genes', 'ExcludedGenes.txt'
)
plot_dir = here('plots', '09_HD_cell_level', 'registration_banksy')

spe = loadHDF5SummarizedExperiment(spe_dir)

#   For both Banksy and FICTURE, read in enrichment modeling results and subset
#   to genes that significantly differentiate a cluster from the others
#   (markers)
banksy_model_results = readRDS(banksy_model_path)$enrichment |>
    as_tibble() |>
    select(ensembl, matches('^fdr_')) |>
    pivot_longer(
        cols = matches('^fdr_'), names_prefix = "fdr_", names_to = 'cluster',
        values_to = 'fdr'
    ) |>
    filter(fdr < 0.05) |>
    mutate(method = "Banksy")

ficture_model_results = readRDS(ficture_model_path)$enrichment |>
    as_tibble() |>
    select(ensembl, matches('^fdr_')) |>
    pivot_longer(
        cols = matches('^fdr_'), names_prefix = "fdr_", names_to = 'cluster',
        values_to = 'fdr'
    ) |>
    filter(fdr < 0.05) |>
    mutate(method = "FICTURE")

#   Read in vector of genes with at least one bad probe and all bad probes,
#   respectively
bad_genes = read_delim(
        bad_path, col_names = c('ensembl', 'symbol'), delim = "|"
    ) |>
    pull(ensembl)
very_bad_genes = read_delim(
        very_bad_path, col_names = c('ensembl', 'symbol'), delim = "|"
    ) |>
    pull(ensembl)

#   Calculate representation of cluster markers in sets of bad genes
model_df = rbind(banksy_model_results, ficture_model_results) |>
    group_by(cluster, method) |>
    slice_head(n = 100) |>
    summarize(
        prop_bad = length(which(ensembl %in% bad_genes)) / n(),
        prop_very_bad = length(which(ensembl %in% very_bad_genes)) / n()
    ) |>
    pivot_longer(
        cols = c(prop_bad, prop_very_bad), names_prefix = 'prop_',
        names_to = 'badness', values_to = 'prop'
    )

#   Calculate y-intercepts representing the proportions expected at random when
#   sampling all genes
int_df = tibble(
    badness = c('bad', 'very_bad'),
    prop = c(length(bad_genes) / nrow(spe), length(very_bad_genes) / nrow(spe))
)

p = ggplot(model_df, aes(x = method, y = prop, color = method)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter() +
    facet_wrap(~badness) +
    geom_hline(data = int_df, aes(yintercept = prop), linetype = 'dashed') +
    coord_cartesian(ylim = c(0, max(model_df$prop))) +
    theme_bw(base_size = 20) +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    ) +
    labs(
        x = 'Gene Set', y = 'Prop. Markers w/ Bad Probes',
        color = "Clustering\nMethod"
    )
pdf(file.path(plot_dir, 'bad_probes.pdf'))
print(p)
dev.off()

session_info()

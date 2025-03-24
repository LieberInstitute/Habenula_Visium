#   10X announced that many probes in the probe set as input to SpaceRanger were
#   "bad". This script investigates the prevalence of bad probes in top markers
#   for Banksy clusters

library(here)
library(tidyverse)
library(spatialLIBD)
library(HDF5Array)
library(sessioninfo)

spe_dir = here('processed-data', '09_HD_cell_level', 'spe_norm_filtered')
model_path = here(
    'processed-data', '09_HD_cell_level', 'registration_banksy',
    'modeling_results', '1.rds'
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

#   Read in enrichment modeling results and subset to genes that significantly
#   differentiate a cluster from the others (markers)
model_results = readRDS(model_path)$enrichment |>
    as_tibble() |>
    select(ensembl, matches('^fdr_')) |>
    pivot_longer(
        cols = matches('^fdr_'), names_prefix = "fdr_", names_to = 'cluster',
        values_to = 'fdr'
    ) |>
    filter(fdr < 0.05)

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
model_df = model_results |>
    filter(fdr < 0.05) |>
    group_by(cluster) |>
    slice_head(n = 100) |>
    summarize(
        prop_bad = length(which(ensembl %in% bad_genes)) / n(),
        prop_very_bad = length(which(ensembl %in% very_bad_genes)) / n()
    )

p = ggplot(model_df, aes(x = 1, y = prop_bad)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter() +
    geom_hline(yintercept = length(bad_genes) /  nrow(spe), linetype = 'dashed') +
    theme_bw(base_size = 20) +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    ) +
    labs(x = '', y = 'Prop. Markers w/ Bad Probes')
pdf(file.path(plot_dir, 'bad_probes.pdf'))
print(p)
dev.off()

session_info()

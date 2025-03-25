#   10X announced that many probes in the probe set as input to SpaceRanger were
#   "bad". This script investigates the prevalence of bad probes in top markers
#   for Banksy and FICTURE clusters

library(here)
library(tidyverse)
library(spatialLIBD)
library(HDF5Array)
library(UpSetR)
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
bs_model_path = here(
    "processed-data", "05_layer_differential_expression",
    "modeling_results_BS", "modeling_results_BayesSpace_k09.Rdata"
)
sn_fine_model_path = here(
    "processed-data", "05_snRNA-seq_model_stats",
    "enrichment_final_Annotations.rds"
)
sn_broad_model_path = here(
    "processed-data", "05_snRNA-seq_model_stats",
    "enrichment_final_Annotations_broad.rds"
)
bad_path = here(
    'processed-data', '10_HD_bin_level', 'bad_probe_genes',
    'OtherGenesWithExcludedProbes.txt'
)
very_bad_path = here(
    'processed-data', '10_HD_bin_level', 'bad_probe_genes', 'ExcludedGenes.txt'
)
banksy_cluster_path = here(
    'processed-data', '09_HD_cell_level', 'banksy', 'lambda0_8',
    'leiden_res1_no_harmony.csv'
)
plot_dir = here('plots', '09_HD_cell_level', 'registration_banksy')

################################################################################
#   Functions
################################################################################

#   Take a dataframe of enrichment modeling results, along with method name,
#   and return a tibble properly formatted for plotting
clean_model_results = function(model_results, method_name) {
    model_results = model_results |>
         as_tibble() |>
        select(ensembl, matches('^fdr_')) |>
        pivot_longer(
            cols = matches('^fdr_'), names_prefix = "fdr_",
            names_to = 'cluster', values_to = 'fdr'
        ) |>
        filter(fdr < 0.05) |>
        mutate(method = method_name)
    
    return(model_results)
}

spe = loadHDF5SummarizedExperiment(spe_dir)

################################################################################
#   Read in and clean enrichment modeling results for each dataset
################################################################################

model_results_list = list()
model_results_list[['Banksy']] = clean_model_results(
    readRDS(banksy_model_path)$enrichment, 'Banksy'
)
model_results_list[['FICTURE']] = clean_model_results(
    readRDS(ficture_model_path)$enrichment, 'FICTURE'
)
model_results_list[['BayesSpace']] = clean_model_results(
    get(load(bs_model_path))$enrichment, 'BayesSpace'
)
model_results_list[['snRNA-seq Fine']] = clean_model_results(
    readRDS(sn_fine_model_path), 'snRNA-seq Fine'
)
model_results_list[['snRNA-seq Broad']] = clean_model_results(
    readRDS(sn_fine_model_path), 'snRNA-seq Broad'
)

################################################################################
#   Read in sets of problematic genes from 10X
################################################################################

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

################################################################################
#   Calculate overlap with bad genes and plot
################################################################################

#   Calculate representation of cluster markers in sets of bad genes
model_df = do.call(rbind, model_results_list) |>
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

#   Box plots
p = ggplot(model_df, aes(x = method, y = prop, color = method)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter() +
    facet_wrap(~badness) +
    geom_hline(data = int_df, aes(yintercept = prop), linetype = 'dashed') +
    coord_cartesian(ylim = c(0, max(model_df$prop))) +
    guides(color = "none") +
    theme_bw(base_size = 20) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(x = 'Clustering Method', y = 'Prop. Markers w/ Bad Probes')
pdf(file.path(plot_dir, 'bad_probes.pdf'), width = 10)
print(p)
dev.off()

#   UpSet plots
gene_sets = lapply(model_results_list, function(x) x$ensembl)
gene_sets[['Bad Genes']] = bad_genes
gene_sets[['Very Bad Genes']] = very_bad_genes
p = upset(
    fromList(gene_sets),
    sets = names(gene_sets),
    order.by = "freq",
    sets.bar.color = "steelblue",
    text.scale = 1,
    number.angles = 15
)
pdf(file.path(plot_dir, 'bad_probes_upset.pdf'), width = 10)
print(p)
dev.off()

#   Plot worst outliers spatially for one sample
cluster_df = read_csv(banksy_cluster_path, show_col_types = FALSE)
stopifnot(identical(spe$key, cluster_df$key))
spe$banksy = case_when(
    cluster_df$banksy_lambda0_8 %in% c(4, 6, 11) ~ 'Habenula',
    cluster_df$banksy_lambda0_8 == 10 ~ 'Top Outlier: Bad',
    cluster_df$banksy_lambda0_8 == 15 ~ 'Top Outlier: Very Bad',
    TRUE ~ 'Other'
)

custom_colors = c(
    'Habenula' = '#96bbbb',
    'Top Outlier: Bad' = '#0125C4',
    'Top Outlier: Very Bad' = '#90121C',
    'Other' = '#F1C606'
)
p = vis_clus(
        spe, sampleid = "H1-MVPY9BW_A1_8433", clustervar = "banksy",
        is_stitched = TRUE, point_size = 20, spatial = FALSE,
        colors = custom_colors
    ) +
    guides(fill = guide_legend(override.aes = list(size = 8)))
png(file.path(plot_dir, 'bad_probes_spatial.png'), width = 1500, height = 1500)
print(p)
dev.off()

session_info()

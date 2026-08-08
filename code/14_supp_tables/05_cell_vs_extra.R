library(tidyverse)
library(here)
library(sessioninfo)
library(qs2)
library(SpatialExperiment)

nuclear_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'contamination',
    'spe', 'raw.qs2'
)
cell_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
extra_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA',
    'extracellular', 'spe_norm_filtered.qs2'
)
plot_path = here('plots', '14_supp_tables', 'cell_vs_extra_boxplots.pdf')

################################################################################
#   Functions
################################################################################

spe_stats = function(spe, dataset_name) {
    tibble(
        sum_umi = colSums(assays(spe)$counts),
        n_genes = colSums(assays(spe)$counts > 0),
        num_bins = spe$bin_count,
        dataset = dataset_name
    )
}

################################################################################
#   Main
################################################################################

spe_nuclear = qs_read(nuclear_path)
spe_cell = readRDS(cell_path)
spe_extra = qs_read(extra_path)

names(assays(spe_nuclear)) = "counts"

#   Subset to the intersection of cells present in all datasets
shared_cells = intersect(
    intersect(colnames(spe_cell), colnames(spe_extra)),
    colnames(spe_nuclear)
)
spe_nuclear = spe_nuclear[, shared_cells]
spe_cell = spe_cell[, shared_cells]
spe_extra = spe_extra[, shared_cells]

metric_names = c(
    "sum_umi" = "Total UMI Count",
    "n_genes" = "Detected Genes",
    "num_bins" = "Num 2um Bins"
)
p = rbind(
        spe_stats(spe_nuclear, "nuclear"),
        spe_stats(spe_cell, "cell"),
        spe_stats(spe_extra, "extra")
    ) |>
    pivot_longer(
        cols = c("sum_umi", "n_genes", "num_bins"),
        names_to = "metric", values_to = "value"
    ) |>
    mutate(
        dataset = factor(dataset, levels = c("nuclear", "cell", "extra")),
        metric = factor(
            metric_names[metric], levels = metric_names
        )
    ) |>
    ggplot(aes(x = dataset, y = value, fill = dataset)) +
        geom_boxplot(outlier.shape = NA) +
        facet_wrap(~metric, scales = "free_y") +
        theme_bw(base_size = 18) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        labs(x = "Dataset", y = "Metric Value") +
        guides(fill = 'none')
pdf(plot_path, height = 5)
print(p)
dev.off()

session_info()

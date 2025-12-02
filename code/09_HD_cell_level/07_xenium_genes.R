#   This script explores expression of the 5001 genes panelled in Xenium in the
#   Visium HD data

library(here)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(sessioninfo)

spe_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'spe_norm_filtered.rds'
)
xenium_path = here(
    'processed-data', '09_HD_cell_level',
    'XeniumPrimeHuman5Kpan_tissue_pathways_metadata.csv'
)
plot_dir = here('plots', '09_HD_cell_level', 'new_samples2', 'QC')

################################################################################
#   Functions
################################################################################

get_exp_genes = function(spe, num_points = 200) {
    stopifnot(length(unique(spe$sample_id)) == 1)

    a = assays(spe)$counts

    #   Ensure we don't have more than one point per integer UMI threshold,
    #   which causes step-like line plots
    num_points = min(max(a) / 3, num_points)

    #   Vectors of UMI cutoffs and proportion of genes with at least one cell
    #   greater than that cutoff, respectively
    cutoff_vals = seq_len(num_points) * max(a) / 3 / num_points
    prop_genes = rep(0, num_points)

    i = 1
    for (umi_cutoff in cutoff_vals) {
        #   For speed, drop genes permanently that don't meet the criterion.
        #   Then count the proportion of the original that are left
        a = a[rowSums(a > umi_cutoff) > 1,,drop = FALSE]
        prop_genes[i] = nrow(a) / 5001
        i = i + 1
    }

    b = tibble(
        umi_cutoff = cutoff_vals,
        prop_genes = prop_genes,
        sample_id = unique(spe$sample_id)
    )

    return(b)
}

################################################################################
#   Main
################################################################################

spe = readRDS(spe_path)

xenium = read_csv(xenium_path)$gene_id
stopifnot(length(xenium) == 5001)

#   Some Xenium genes simply aren't measured in Visium HD
num_missing = sum(!(xenium %in% rownames(spe)))
message(
    sprintf(
        "%s of %s Xenium genes are not measured in Visium HD (%s%%)",
        num_missing,
        length(xenium),
        round(100 * num_missing / length(xenium), 1)
    )
)

#   Subset to Xenium genes
spe = spe[xenium[xenium %in% rownames(spe)],]

exp_df = do.call(
    rbind, 
    lapply(
        unique(spe$sample_id),
        function(x) {
            get_exp_genes(spe[, spe$sample_id == x])
        }
    )
)

#   Line plot colored by region of proportion of genes having at least one cell
#   with counts greater than various thresholds
p = exp_df |>
    filter(prop_genes > 0.01) |>
    ggplot(
            aes(
                x = umi_cutoff, y = prop_genes, color = sample_id,
                group = sample_id
            )
        ) +
        geom_line() +
        theme_bw(base_size = 25) +
        labs(x = "UMI cutoff", y = "Prop. of genes", color = "Sample ID")

pdf(file.path(plot_dir, 'xenium_exp_genes_umi_cutoffs.pdf'), width = 10)
print(p)
dev.off()

session_info()

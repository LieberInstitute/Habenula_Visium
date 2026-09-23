#   High-resolution spot plots of genes profiled by RNAscope. Will need this
#   main figure 1

library(here)
library(tidyverse)
library(spatialLIBD)
library(sessioninfo)

sample_id = 'Br9090_1'
spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'no_secondary', 'misc_paper_figs'
)
all_genes = c("TAC3", "GPR151", "POU4F1", "MBP")
num_pix = 3000

spe = readRDS(spe_path)

all_genes_ensembl = rowData(spe)$gene_id[
    match(all_genes, rowData(spe)$gene_name)
]
stopifnot(!any(is.na(all_genes_ensembl)))

for (i in seq_along(all_genes)) {
    #   Run twice to overcome a bug with different behavior on the first plot
    for (j in seq_len(2)) {
        p = vis_gene(
                spe, sampleid = sample_id, geneid = all_genes_ensembl[i],
                spatial = FALSE, is_stitched = TRUE, point_size = 50,
                cap_percentile = 0.995
            ) +
            theme_void(base_size = 50) +
            theme(
                legend.key.size = unit(2, "cm"), legend.direction = "horizontal"
            )
    }
  
    png(
        file.path(plot_dir, sprintf('%s.png', all_genes[i])),
        width = num_pix, height = num_pix
    )
    print(p)
    dev.off()
}

session_info()

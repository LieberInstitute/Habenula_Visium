#   Plot top markers for various interneuron subtypes

library(here)
library(tidyverse)
library(spatialLIBD)
library(sessioninfo)
library(readxl)
library(cowplot)

spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
marker_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'supp_table_1_filt.xlsx'
)
plot_dir = here('plots', '09_HD_cell_level', 'no_secondary', 'interneurons')
top_n = 10
sample_ids = c('Br9090_1', 'Br8433_1')

dir.create(plot_dir, showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

vis_gene_clean = function(spe, gene_vec, sample_id) {
    #   Run twice to overcome a bug with different behavior on the first plot
    for (i in seq_len(2)) {
        p = vis_gene(
            spe, sampleid = sample_id, geneid = gene_vec,
            is_stitched = TRUE, point_size = 20, spatial = FALSE,
            cap_percentile = 0.99
        )
    }
    return(p)
}

################################################################################
#   Main
################################################################################

spe = readRDS(spe_path)

marker_df = read_excel(marker_path) |>
    mutate(
        gene_id = rowData(spe)$gene_id[match(gene, rowData(spe)$gene_name)]
    ) |>
    filter(!is.na(gene_id), fold_change > 1) |>
    group_by(gene_id) |>
    filter(n() == 1) |>
    ungroup() |>
    group_by(cell_type) |>
    arrange(desc(auroc)) |>
    slice_head(n = top_n) |>
    ungroup() |>
    select(gene_id, cell_type)

for (cell_type in unique(marker_df$cell_type)) {
    gene_vec = marker_df |>
        filter(cell_type == !!cell_type) |>
        pull(gene_id)

    plot_list = list()
    for (sample_id in sample_ids) {
        plot_list[[sample_id]] = vis_gene_clean(
            spe, gene_vec = gene_vec, sample_id = sample_id
        )
    }

    png(
        file.path(plot_dir, sprintf("%s.png", cell_type)),
        width = 3000, height = 1500
    )
    print(plot_grid(plotlist = plot_list), nrow = 1)
    dev.off()
}

session_info()

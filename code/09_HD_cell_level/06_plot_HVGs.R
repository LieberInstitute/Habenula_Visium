library(here)
library(tidyverse)
library(SpatialExperiment)
library(HDF5Array)
library(spatialLIBD)
library(sessioninfo)
library(cowplot)
library(viridis)

spe_dir = here('processed-data', '09_HD_cell_level', 'spe_norm_filtered')
plot_dir = here('plots', '09_HD_cell_level', 'HVGs')
hvg_path = here('processed-data', '09_HD_cell_level', 'HVGs.txt')
top_n = 12
px_per_plot = 300

dir.create(plot_dir, showWarnings = FALSE)

plot_one_sample = function(sample_id, spe, hvg) {
    #   Subset to this sample
    spe = spe[, spe$sample_id == sample_id]

    #   Cut off expression at the top 5% to make the color range more dynamic
    hvg_exp = as.matrix(assays(spe)$logcounts[hvg,])
    for (this_gene in hvg) {
        cutoff = sort(hvg_exp[this_gene,], decreasing = TRUE)[
            as.integer(5 * ncol(spe) / 100)
        ]
        assays(spe)$logcounts[this_gene,] = pmin(hvg_exp[this_gene,], cutoff)
    }

    #   Note: a lot of manual plotting code is used in place of vis_gene because
    #   cropping the in-tissue spots for this cell-level Visium HD data is not
    #   currently possible through vis_gene, leading to excessive whitespace

    #   Gather expression and spatial coordinates into a tidy tibble
    exp_df = assays(spe)$logcounts[hvg, ] |>
        t() |>
        as_tibble() |>
        cbind(spatialCoords(spe)) |>
        as_tibble() |>
        mutate(
            x = pxl_col_in_fullres,
            y = max(pxl_row_in_fullres) - pxl_row_in_fullres
        )

    plot_list = list()
    for (i in seq_len(top_n)) {
        plot_list[[i]] = ggplot(
                exp_df,
                aes(x = x, y = y, color = !!sym(hvg[i]), fill = !!sym(hvg[i]))
            ) +
            geom_point(size = 0.0003) +
            scale_fill_viridis() +
            scale_color_viridis() +
            labs(
                fill = 'expr',
                title = rowData(spe)$gene_name[match(hvg[i], rownames(spe))]
            ) +
            guides(color = "none") +
            theme_bw(base_size = 15) +
            #   Remove pretty much everything related to x- and y-axis labels
            theme(
                axis.title.x = element_blank(), axis.title.y = element_blank(),
                axis.text.x = element_blank(), axis.text.y = element_blank(),
                axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
                plot.title = element_text(size = 25)
            )
    }

    png(
        file.path(plot_dir, sprintf('top_%s_HVGs_%s.png', top_n, sample_id)),
        width = px_per_plot * top_n / 3, height = px_per_plot * 3
    )
    print(plot_grid(plotlist = plot_list, nrow = 3))
    dev.off()
}

hvg = readLines(hvg_path)[seq_len(top_n)]
spe = loadHDF5SummarizedExperiment(spe_dir)
spe$exclude_overlapping = FALSE

for (sample_id in unique(spe$sample_id)) {
    message("Plotting HVGs for sample ", sample_id)
    plot_one_sample(sample_id, spe, hvg)
}

session_info()

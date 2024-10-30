library(here)
library(tidyverse)
library(SpatialExperiment)
library(HDF5Array)
library(spatialLIBD)
library(sessioninfo)
library(cowplot)
library(viridis)

spe_dir = here('processed-data', '10_HD_bin_level', 'spe_norm')
plot_dir = here('plots', '10_HD_bin_level')
svg_path = here(
    'processed-data', '10_HD_bin_level', 'nnSVG_out', 'H1-W369TJK_D1_9090.csv'
)
top_n = 12
px_per_plot = 500

svg = read_csv(svg_path, show_col_types = FALSE) |>
    arrange(rank) |>
    slice_head(n = top_n) |>
    pull(gene_id)

spe = loadHDF5SummarizedExperiment(spe_dir)
spe$exclude_overlapping = FALSE

stopifnot(all(svg %in% rownames(spe)))

#   Note: a lot of manual plotting code is used in place of vis_gene because
#   cropping the in-tissue spots for this Visium HD data is not
#   currently possible through vis_gene, leading to excessive whitespace

#   Gather expression and spatial coordinates into a tidy tibble
exp_df = assays(spe)$logcounts[svg, ] |>
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
            aes(x = x, y = y, color = !!sym(svg[i]), fill = !!sym(svg[i]))
        ) +
        geom_point(size = 0.0003) +
        scale_fill_viridis() +
        scale_color_viridis() +
        labs(
            fill = 'expr',
            title = rowData(spe)$gene_name[match(svg[i], rownames(spe))]
        ) +
        guides(color = "none") +
        #   Remove pretty much everything related to x- and y-axis labels
        theme(
            axis.title.x = element_blank(), axis.title.y = element_blank(),
            axis.text.x = element_blank(), axis.text.y = element_blank(),
            axis.ticks.x = element_blank(), axis.ticks.y = element_blank()
        )
}

png(
    file.path(plot_dir, sprintf('top_%s_nnSVGs.png', top_n)),
    width = px_per_plot * top_n / 3, height = px_per_plot * 3
)
plot_grid(plotlist = plot_list, nrow = 3)
dev.off()

session_info()

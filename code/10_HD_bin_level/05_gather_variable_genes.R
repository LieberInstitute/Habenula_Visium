library(here)
library(tidyverse)
library(SpatialExperiment)
library(HDF5Array)
library(spatialLIBD)
library(sessioninfo)
library(cowplot)
library(viridis)
library(scran)

spe_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
plot_dir = here('plots', '10_HD_bin_level', 'no_secondary', 'variable_genes')
svg_paths = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'nnSVG_out', '%s.csv'
)
svg_path_out = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'nnSVG_out',
    'merged_SVGs.txt'
)
hvg_path = here('processed-data', '09_HD_cell_level', 'new_samples2', 'HVGs.txt')

num_svg = 1000
top_n = 12
px_per_plot = 300

dir.create(plot_dir, showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

plot_VGs = function(spe, genes, plot_str, percentile = 0.98) {
    #   Cutoff expression at the top to make the color range more dynamic
    vg_exp = as.matrix(assays(spe)$logcounts[genes,])
    for (this_gene in genes) {
        cutoff = quantile(vg_exp[this_gene,], percentile)
        cutoff = ifelse(cutoff == 0, max(vg_exp[this_gene,]), cutoff)
        
        vg_exp[this_gene,] = pmin(vg_exp[this_gene,], cutoff)
    }

    #   Note: a lot of manual plotting code is used in place of vis_gene because
    #   cropping the in-tissue spots for this Visium HD data is not
    #   currently possible through vis_gene, leading to excessive whitespace

    #   Gather expression and spatial coordinates into a tidy tibble
    exp_df = vg_exp |>
        t() |>
        as_tibble() |>
        cbind(spatialCoords(spe)) |>
        as_tibble() |>
        mutate(
            x = pxl_col_in_fullres,
            y = max(pxl_row_in_fullres) - pxl_row_in_fullres
        )

    plot_list = list()
    for (i in seq_len(length(genes))) {
        plot_list[[i]] = ggplot(
                exp_df,
                aes(
                    x = x, y = y, color = !!sym(genes[i]),
                    fill = !!sym(genes[i])
                )
            ) +
            geom_point(size = 0.1) +
            scale_fill_viridis() +
            scale_color_viridis() +
            labs(
                fill = 'expr',
                title = rowData(spe)$gene_name[match(genes[i], rownames(spe))]
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
        file.path(plot_dir, sprintf('top_%s_%s.png', top_n, plot_str)),
        width = px_per_plot * top_n / 3, height = px_per_plot * 3
    )
    print(plot_grid(plotlist = plot_list, nrow = 3))
    dev.off()
}

################################################################################
#   Read in HVGs and calculate SVGs
################################################################################

spe = readRDS(spe_path)
spe$exclude_overlapping = FALSE

hvg = readLines(hvg_path)

#   Read in SVGs for each sample
svg_list = list()
for (sample_id in unique(spe$sample_id)) {
    svg_list[[sample_id]] = sprintf(svg_paths, sample_id) |>
        read_csv(show_col_types = FALSE) |>
        select(gene_id, rank) |>
        mutate(sample_id = sample_id)
}

#   Quite a few genes were considered (passed expression cutoffs) in all
#   samples, so the "filter(n() == 5)" step is not merely selecting for
#   high-expression genes
svg = do.call(rbind, svg_list) |>
    group_by(gene_id) |>
    filter(n() == 5)

message(
    sprintf(
        "%s unique genes were considered as candidate SVGs in all 5 samples",
        svg |>
            pull(gene_id) |>
            unique() |>
            length()        
    )
)

#   Take top SVGs by average rank across samples
svg = svg |>
    summarize(avg_rank = mean(rank)) |>
    arrange(avg_rank) |>
    slice_head(n = num_svg) |>
    pull(gene_id)

writeLines(svg, svg_path_out)

################################################################################
#   Plot HVGs and SVGs
################################################################################

for (sample_id in unique(spe$sample_id)) {
    small_spe = spe[, spe$sample_id == sample_id]
    plot_VGs(small_spe, hvg[seq_len(top_n)], sprintf("HVGs_%s", sample_id))
    plot_VGs(small_spe, svg[seq_len(top_n)], sprintf("SVGs_%s", sample_id))
}

session_info()

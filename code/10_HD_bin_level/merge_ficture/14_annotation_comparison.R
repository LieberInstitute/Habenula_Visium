library(here)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(HDF5Array)
library(sessioninfo)
library(cowplot)

spe_dir = here('processed-data', '09_HD_cell_level', 'spe_norm')
cor_path = here(
    'processed-data', '09_HD_cell_level', 'banksy', 'cor_vs_snRNA-seq.rds'
)
singler_broad_path = here(
    'processed-data', '09_HD_cell_level', 'singler', 'broad.csv'
)
singler_fine_path = here(
    'processed-data', '09_HD_cell_level', 'singler', 'fine.csv'
)
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'banksy', 'k%s.csv'
)
plot_dir = here('plots', '09_HD_cell_level', 'banksy')
px_per_plot = 500

parse_cor_mat = function(cor_list, res) {
    cor_df = do.call(rbind, cor_list[[res]]) |>
        as.data.frame()
    cor_df$cell_type = apply(
            cor_df, 1, function(x) { colnames(cor_df)[which.max(x)] }
        ) |>
        unname()
    cor_df = cor_df |>
        rownames_to_column('rn') |>
        mutate(
            k = str_extract(rn, '^k([0-9]+)_', group = 1) |>
                as.numeric(),
            cluster_num = str_extract(rn, '_([0-9]+) ', group = 1) |>
                as.numeric(),
            res = {{ res }}
        ) |>
        as_tibble() |>
        select(k, res, cluster_num, cell_type)
    
    return(cor_df)
}

#   This script will include two plots:
#   spot plot pair for k = 8 broad (singleR vs spatial reg)
#   line plot (x = k, y = agreement, color = resolution)

spe = loadHDF5SummarizedExperiment(spe_dir)

#   Add cell types called by SingleR to the SPE
spe$singler_broad = read_csv(singler_broad_path, show_col_types = FALSE) |>
    pull(labels) |>
    factor()
spe$singler_fine = read_csv(singler_fine_path, show_col_types = FALSE) |>
    pull(labels) |>
    factor()

#   Read banksy clustering results into the SPE
for (k in 2:28) {
    cluster_df = read_csv(sprintf(cluster_path, k), show_col_types = FALSE)

    stopifnot(setequal(as.numeric(colnames(spe)), cluster_df$key))
    spe[[paste0('banksy_k', k)]] = cluster_df$banksy_lambda0.2[
        match(as.numeric(colnames(spe)), cluster_df$key)
    ]
}

#   Associate a single cell type with each cluster value for each k at both
#   resolutions
cor_list = readRDS(cor_path)
cor_df = rbind(
    parse_cor_mat(cor_list, 'broad'), parse_cor_mat(cor_list, 'fine')
)

#   Annotate banksy cluster results with cell types in the SPE
conc_df_list = list()
for (res in c('broad', 'fine')) {
    for (k in 2:28) {
        this_cor_df = cor_df |>
            filter(k == {{ k }}, res == {{ res }})
        
        spe[[sprintf('anno_k%s_%s', k, res)]] = this_cor_df$cell_type[
            match(spe[[paste0('banksy_k', k)]], this_cor_df$cluster_num)
        ]

        conc_df_list[[paste0(k, res)]] = tibble(
            k = k,
            res = res,
            concordance = mean(
                spe[[sprintf('anno_k%s_%s', k, res)]] == spe[[paste0('singler_', res)]]
            )
        )
    }
}
conc_df = do.call(rbind, conc_df_list)

#   Explore agreement of spatial registration with SingleR results at different
#   k values and cell-type resolutions
p = ggplot(conc_df, aes(x = k, y = concordance, color = res, group = res)) +
    geom_line() +
    scale_y_continuous(limits = c(0, max(conc_df$concordance))) +
    scale_x_continuous(breaks = seq_len(14) * 2) +
    theme_bw(base_size = 15) +
    labs(x = 'Banksy k value', y = '% agreement', color = 'Cell-type\nresolution')
pdf(file.path(plot_dir, 'annotation_concordance.pdf'), width = 8, height = 6)
print(p)
dev.off()

#   Define a named vector of colors for broad and fine cell types
cell_colors = list()

temp = cor_df |>
    filter(res == 'broad') |>
    pull(cell_type) |>
    unique() |>
    sort()
cell_colors[['broad']] = rainbow(length(temp))
names(cell_colors[['broad']]) = temp

temp = cor_df |>
    filter(res == 'fine') |>
    pull(cell_type) |>
    unique() |>
    sort()
cell_colors[['fine']] = rainbow(length(temp))
names(cell_colors[['fine']]) = temp

#   Note: a lot of manual plotting code is used in place of vis_clus because
#   cropping the in-tissue spots for this cell-level Visium HD data is not
#   currently possible through vis_clus, leading to excessive whitespace

#   Gather cluster results and spatial coordinates into a tidy tibble
col_names = c(
    sprintf('anno_k%s_broad', length(cell_colors[['broad']])),
    sprintf('anno_k%s_fine', length(cell_colors[['fine']])),
    'singler_broad',
    'singler_fine'
)
exp_df = colData(spe)[, col_names] |>
    cbind(spatialCoords(spe)) |>
    as_tibble() |>
    mutate(
        x = pxl_col_in_fullres,
        y = max(pxl_row_in_fullres) - pxl_row_in_fullres
    )

for (res in c('broad', 'fine')) {
    plot_list = list()
    for (col_name in col_names[grep(res, col_names)]) {
        plot_list[[col_name]] = ggplot(
                exp_df,
                aes(
                    x = x, y = y, color = !!sym(col_name),
                    fill = !!sym(col_name)
                )
            ) +
            geom_point(size = 0.001) +
            scale_color_manual(values = cell_colors[[res]]) +
            labs(color = 'cell type', title = col_name) +
            guides(
                color = guide_legend(override.aes = list(size = 4)),
                fill = "none"
            ) +
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
        file.path(plot_dir, sprintf('annotation_spot_plots_%s.png', res)),
        width = px_per_plot * 2, height = px_per_plot
    )
    print(plot_grid(plotlist = plot_list, nrow = 1))
    dev.off()
}

session_info()
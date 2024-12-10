library(here)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(HDF5Array)
library(sessioninfo)

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
        
        spe[[sprintf('anno_k%s_%s', k, res)]] = cor_df$cell_type[
            match(spe[[paste0('banksy_k', k)]], cor_df$cluster_num)
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

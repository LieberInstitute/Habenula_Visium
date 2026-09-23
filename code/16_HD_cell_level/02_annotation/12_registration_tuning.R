#   Spatial registration results seem sensitive to the top_n parameter in
#   layer_stat_cor. Here we try a range of values for this parameter, rating
#   the quality of the results by how Banksy-multiome correlation heatmaps and
#   Banksy-Yalcinbas heatmaps correlate (higher is better, since for the most
#   part the cell-type definitions match).

library(here)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(sessioninfo)

all_resolutions = c(seq_len(20) / 10, 4, 8)
plot_path = here(
    'plots', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'tuning_top_n.pdf'
)
model_paths = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'modeling_results',
    sprintf('%s.rds', c(sub('\\.', '_', as.character(all_resolutions))))
)
ref_yalcinbas_path = here(
    "processed-data", "05_snRNA-seq_model_stats",
    "enrichment_final_Annotations.rds"
)
ref_multiome_path = '/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/11_link_prep/04_registration_wrapper/model_results.rds'
all_top_n = c(10, 30, 100, 300, 1000, 3000, 10000, 13000, Inf)
multiome_clusters = c(
    "Astrocyte", "Endo", "Excit.Thal", "LHb.1.3.4", "LHb.2.7", "MHb.1", "MHb.3",
    "Microglia", "Oligo", "OPC"
)

dir.create(dirname(plot_path), showWarnings = FALSE)

#   Read in enrichment stats for Banksy clusters at all Leiden resolutions
hd_stats = lapply(model_paths, function(path) readRDS(path)$enrichment)
names(hd_stats) = as.character(all_resolutions)

multiome_stats = list(
    enrichment = readRDS(ref_multiome_path)$enrichment |>
        filter(!duplicated(ensembl))
)
yalcinbas_stats = list(enrichment = readRDS(ref_yalcinbas_path))

#   Genes can be filtered, leading to different numbers for each resolution,
#   but for plotting purposes we take the mean to get a general idea
num_genes_total = mean(sapply(hd_stats, function(x) length(unique(x$ensembl))))

cor_df_list = list()
for (this_top_n in all_top_n) {
    for (this_resolution in all_resolutions) {
        if (this_top_n == Inf) {
            multiome_cor = layer_stat_cor(
                hd_stats[[as.character(this_resolution)]],
                modeling_results = multiome_stats, model_type = "enrichment"
            )
            yalcinbas_cor = layer_stat_cor(
                hd_stats[[as.character(this_resolution)]],
                modeling_results = yalcinbas_stats, model_type = "enrichment"
            )
        } else {
            multiome_cor = layer_stat_cor(
                hd_stats[[as.character(this_resolution)]],
                modeling_results = multiome_stats, model_type = "enrichment",
                top_n = this_top_n
            )
            yalcinbas_cor = layer_stat_cor(
                hd_stats[[as.character(this_resolution)]],
                modeling_results = yalcinbas_stats, model_type = "enrichment",
                top_n = this_top_n
            )
        }

        multiome_df = multiome_cor |>
            as.data.frame() |>
            rownames_to_column('hd_cluster') |>
            as_tibble() |>
            pivot_longer(
                -hd_cluster, names_to = 'ref_cluster',
                values_to = 'multiome_cor_val'
            ) |>
            filter(ref_cluster %in% multiome_clusters) |>
            mutate(resolution = this_resolution, top_n = this_top_n)
        yalcinbas_df = yalcinbas_cor |>
            as.data.frame() |>
            rownames_to_column('hd_cluster') |>
            as_tibble() |>
            pivot_longer(
                -hd_cluster, names_to = 'ref_cluster',
                values_to = 'yalcinbas_cor_val'
            ) |>
            mutate(
                ref_cluster = case_when(
                    grepl('^LHb\\.[134]', ref_cluster) ~ 'LHb.1.3.4',
                    grepl('^LHb\\.[27]', ref_cluster) ~ 'LHb.2.7',
                    TRUE ~ ref_cluster
                )
            ) |>
            group_by(hd_cluster, ref_cluster) |>
            summarize(yalcinbas_cor_val = mean(yalcinbas_cor_val)) |>
            ungroup() |>
            filter(ref_cluster %in% multiome_clusters) |>
            mutate(resolution = this_resolution, top_n = this_top_n)
        stopifnot(setequal(multiome_df$ref_cluster, yalcinbas_df$ref_cluster))

        cor_df_list[[length(cor_df_list) + 1]] = inner_join(
                multiome_df, yalcinbas_df,
                by = c("hd_cluster", "ref_cluster", "resolution", "top_n")
            )
    }
}

p = bind_rows(cor_df_list) |>
    group_by(resolution, top_n) |>
    summarize(this_cor = cor(multiome_cor_val, yalcinbas_cor_val)) |>
    group_by(top_n) |>
    summarize(this_cor = mean(this_cor)) |>
    mutate(top_n = ifelse(top_n == Inf, num_genes_total, top_n)) |>
    ggplot(aes(x = top_n, y = this_cor)) +
    geom_point() +
    geom_line() +
    scale_x_log10() +
    labs(
        x = "Number of genes used for correlation",
        y = "Correlation between heatmaps"
    ) +
    theme_bw(base_size = 20)
pdf(plot_path)
print(p)
dev.off()

session_info()

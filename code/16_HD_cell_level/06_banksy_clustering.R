library(here)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(sessioninfo)
library(Banksy)
library(cowplot)
library(scater)

lambda = 0.2
res = c(seq_len(20) / 10, 4, 8)[as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))]
res_neat = paste0('res', sub('\\.', '_', as.character(res)))

spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'banksy',
    'spe_banksy.rds'
)
spe_orig_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
out_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'banksy',
    sprintf('leiden_%s.csv', res_neat)
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'no_secondary', 'banksy',
    sprintf('leiden_%s', res_neat)
)
random_seed = 0

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

spe = readRDS(spe_path)

message(Sys.time(), ' | Performing clustering')
spe = clusterBanksy(
    spe, use_agf = TRUE, lambda = lambda, seed = random_seed,
    algo = "leiden", resolution = res, dimred = "HARMONY"
)

#   Get the names of the cluster and UMAP columns from a more general regex
cluster_name = colnames(colData(spe))[
    grep('^clust_HARMONY', colnames(colData(spe)))
]
rd_name = reducedDimNames(spe)[
    grep(sprintf('^UMAP.*lam%s', lambda), reducedDimNames(spe))
]

#   Fix spatial coordinates (unstagger) for plotting
spe_orig = readRDS(spe_orig_path)
spatialCoords(spe) = spatialCoords(spe_orig[, spe$key])
rm(spe_orig)
gc()

#   Plot clusters and colored UMAP for each sample
for (sample_id in unique(spe$sample_id)) {
    if (length(unique(spe[[cluster_name]])) <= 36) {
        #   First plot the clusters spatially
        p = vis_clus(
                spe, sampleid = sample_id, clustervar = cluster_name,
                is_stitched = TRUE, point_size = 20, spatial = FALSE
            ) +
            guides(fill = guide_legend(override.aes = list(size = 8)))
        png(
            file.path(plot_dir, sprintf('clusters_%s.png', sample_id)),
            width = 1500, height = 1500
        )
        print(p)
        dev.off()
    }

    #   Then UMAP colored by cluster
    p = plotReducedDim(
            spe[,spe$sample_id == sample_id], dimred = rd_name,
            point_size = 0.6, colour_by = cluster_name
        ) +
        theme_bw(base_size = 20) +
        guides(color = guide_legend(override.aes = list(size = 4)))
    png(
        file.path(plot_dir, sprintf('UMAP_%s.png', sample_id)),
        width = 1500, height = 1500
    )
    print(p)
    dev.off()
}

#   Export clusters to CSV
cluster_df = colData(spe) |>
    as_tibble() |>
    select(key, sym(cluster_name))
colnames(cluster_df) = c('key', 'banksy')

write_csv(cluster_df, out_path)

session_info()

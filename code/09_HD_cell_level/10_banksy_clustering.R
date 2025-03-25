library(here)
library(tidyverse)
library(SpatialExperiment)
library(HDF5Array)
library(spatialLIBD)
library(sessioninfo)
library(Banksy)
library(cowplot)
library(scater)
library(getopt)

# Import command-line parameters
spec <- matrix(
    c(
        c("res", "lambda"),
        c("r", "l"),
        rep("1", 2),
        rep("numeric", 2),
        c('Resolution for leiden clustering', 'Banksy hyperparameter lambda')
    ),
    ncol = 5
)
opt <- getopt(spec)

message("Using the following parameters:")
print(opt)

lambda_neat = paste0('lambda', sub('\\.', '_', as.character(opt$lambda)))
res_neat = paste0('res', sub('\\.', '_', as.character(opt$res)))

spe_dir = here(
    'processed-data', '09_HD_cell_level', sprintf('spe_banksy_%s', lambda_neat)
)
out_path = here(
    'processed-data', '09_HD_cell_level', 'banksy', lambda_neat,
    sprintf('leiden_%s.csv', res_neat)
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'banksy', lambda_neat,
    sprintf('leiden_%s', res_neat)
)
random_seed = 0

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

spe = loadHDF5SummarizedExperiment(spe_dir)

message(Sys.time(), ' | Performing clustering')
spe = clusterBanksy(
    spe, use_agf = TRUE, lambda = opt$lambda, seed = random_seed,
    algo = "leiden", resolution = opt$res, dimred = "HARMONY"
)

#   Get the names of the cluster and UMAP columns from a more general regex
cluster_name = colnames(colData(spe))[
    grep(sprintf('^clust.*lam%s', opt$lambda), colnames(colData(spe)))
]
rd_name = reducedDimNames(spe)[
    grep(sprintf('^UMAP.*lam%s', opt$lambda), reducedDimNames(spe))
]

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
colnames(cluster_df) = c('key', sprintf('banksy_%s', lambda_neat))

write_csv(cluster_df, out_path)

session_info()

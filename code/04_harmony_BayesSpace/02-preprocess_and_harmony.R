## Required libraries
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(tidyverse)
library(sessioninfo)
library(scran) ## requires uwot for UMAP
library(uwot)
library(scater)
library(BiocParallel)
library(parallel)
library(PCAtools)
library(ggplot2)
library(Polychrome)
library(harmony)
library(HDF5Array)
library(spatialNAcUtils)


# Define harmony function that would allow us to specify which reduction to use with SingleCellExperiment object
RunHarmony_mod <- function(
    object,
    group.by.vars,
    reduction.use = "PCA",
    dims.use = NULL,
    verbose = TRUE,
    reduction.save = "HARMONY",
    ...
) {

    ## Get PCA embeddings
    if (!"PCA" %in% SingleCellExperiment::reducedDimNames(object)) {
        stop("PCA must be computed before running Harmony.")
    }
    pca_embedding <- SingleCellExperiment::reducedDim(object, reduction.use)
    if (is.null(dims.use)) {
        dims.use <- seq_len(ncol(pca_embedding))
    }

    if (is.null(dims.use)) {
        dims.use <- seq_len(ncol(pca_embedding))
    }
    dims_avail <- seq_len(ncol(pca_embedding))
    if (!all(dims.use %in% dims_avail)) {
        stop("trying to use more dimensions than computed with PCA. Rerun
            PCA with more dimensions or use fewer PCs")
    }

    metavars_df <- SingleCellExperiment::colData(object)
    if (!all(group.by.vars %in% colnames(metavars_df))) {
        stop('Trying to integrate over variables missing in colData')
    }

    harmonyEmbed <- RunHarmony(
        data_mat = pca_embedding[, dims.use],
        meta_data = metavars_df,
        vars_use = group.by.vars,
        return_object = FALSE,
        verbose = verbose,
        ...
    )


    rownames(harmonyEmbed) <- row.names(metavars_df)
    colnames(harmonyEmbed) <- paste0(reduction.save, "_", seq_len(ncol(harmonyEmbed)))
    SingleCellExperiment::reducedDim(object, reduction.save) <- harmonyEmbed

    return(object)
}

tsne_perplex_vals = c("05", "20", "50", "80")
num_cores = detectCores() - 1

## Create output directories
dir_plots <- here("plots", "05_harmony_BayesSpace")
dir_rdata <- here("processed-data", "05_harmony_BayesSpace")
filtered_hdf5_dir = here(
    'processed-data', '05_harmony_BayesSpace', 'spe_filtered_hdf5'
)
harmony_hdf5_dir = here(
    'processed-data', '05_harmony_BayesSpace', 'spe_harmony'
)

dir.create(dir_plots, showWarnings = FALSE)
dir.create(dir_rdata, showWarnings = FALSE)

dir.create(file.path(dir_rdata, "clusters_graphbased"), showWarnings = FALSE)
dir.create(file.path(dir_rdata, "clusters_graphbased_cut_at"), showWarnings = FALSE)

set.seed(20230712)

## Load the data
spe = loadHDF5SummarizedExperiment(filtered_hdf5_dir)

## Plot initial low-dimensional representations prior to batch correction
plotReducedDim(spe, dimred = "PCA", ncomponents = 3, colour_by = "donor")
plotReducedDim(spe, dimred = "GLMPCA_approx", ncomponents = 3, colour_by = "donor")
plotReducedDim(spe, dimred = "PCA", ncomponents = 3, colour_by = "scran_discard")
plotReducedDim(spe, dimred = "GLMPCA_approx", ncomponents = 3, colour_by = "scran_discard")

ggcells(spe, aes(x = GLMPCA_approx.1, y = GLMPCA_approx.2, colour = scran_discard)) +
  geom_point(size = 0.5) +
  facet_wrap(~ slide_num) +
  labs(x = "GLMPC1", y = "GLMPC2", colour = "Discard") + theme_classic()

ggcells(spe, aes(x = GLMPCA_approx.1, y = GLMPCA_approx.2, colour = scran_discard)) +
  geom_point(size = 0.5) +
  facet_wrap(~ sample_id_original) +
  labs(x = "GLMPC1", y = "GLMPC2", colour = "Discard") + theme_classic()



## Perform harmony batch correction
message("Running RunHarmony()")
Sys.time()
spe <- RunHarmony_mod(spe, group.by.vars = c("donor", "slide_num"), verbose = TRUE, plot_convergence = TRUE, kmeans_init_nstart=100, kmeans_init_iter_max=1000)

spe <- RunHarmony_mod(spe, "donor", verbose = TRUE, plot_convergence = TRUE, reduction.use = "PCA", reduction.save = "harmony_donor", kmeans_init_nstart=100, kmeans_init_iter_max=1000)
spe <- RunHarmony_mod(spe, group.by.vars = c("donor", "slide_num"), verbose = TRUE, reduction.save = "harmony_donor_slide")
#spe <- RunHarmony_mod(spe, group.by.vars = c("donor", "sample_id_original"), verbose = TRUE, reduction.save = "HARMONY", plot_convergence = TRUE, lambda = NULL, max_iter = 30)
Sys.time()

#   Perform dimensionality reduction using both PCA and harmony's reduced
#   dimensions
for (dimred_var in c("PCA", "HARMONY")) {
    #   Run TSNE with several perplexity values
    for (perplex in tsne_perplex_vals) {
        message(
            sprintf(
                "Running runTSNE() perplexity %s on %s dimensions",
                perplex, dimred_var
            )
        )
        Sys.time()
        spe <-
            runTSNE(
                spe,
                dimred = dimred_var,
                name = sprintf("TSNE_perplexity%s.%s", perplex, dimred_var),
                perplexity = as.integer(perplex)
            )
        Sys.time()
    }

    #   Explore TSNE results via plots
    pdf(
        file = file.path(
            dir_plots,
            sprintf("tSNE_perplexity%s_%s_sample_id.pdf", perplex, dimred_var)
        ),
        width = 9
    )
    ggplot(
            data.frame(
                reducedDim(
                    spe, sprintf("TSNE_perplexity%s.%s", perplex, dimred_var)
                )
            ),
            aes(x = TSNE1, y = TSNE2, color = factor(spe$sample_id))
        ) +
        geom_point() +
        labs(color = "sample_id") +
        theme_bw()
    dev.off()

    #   Also run UMAP
    message(sprintf("Running runUMAP() on %s dimensions", dimred_var))
    Sys.time()
    spe <- runUMAP(
        spe, dimred = dimred_var, name = sprintf("UMAP.%s", dimred_var),
        BPPARAM = MulticoreParam(num_cores)
    )
    Sys.time()

    #   Explore UMAP results, coloring by both donor and sample ID
    for (color_var in c("sample_id", "donor")) {
        pdf(
            file = file.path(
                dir_plots, sprintf("UMAP_%s_%s.pdf", color_var, dimred_var)
            )
        )
        ggplot(
            data.frame(reducedDim(spe, sprintf("UMAP.%s", dimred_var))),
            aes(x = UMAP1, y = UMAP2, color = factor(spe[[color_var]]))
        ) +
            geom_point() +
            labs(color = color_var) +
            theme_bw()
        dev.off()
    }
}

## Perform graph-based clustering on batch corrected-data
message("Running buildSNNGraph() on HARMONY dimensions")
Sys.time()
g_k10 <- buildSNNGraph(spe, k = 10, use.dimred = "HARMONY")
Sys.time()
save(g_k10, file = file.path(dir_rdata, "g_k10_harmony.Rdata"))

message("Running cluster_walktrap()")
Sys.time()
g_walk_k10 <- igraph::cluster_walktrap(g_k10)
Sys.time()
save(g_walk_k10, file = file.path(dir_rdata, "g_walk_k10_harmony.Rdata"))

clust_k10 <- sort_clusters(g_walk_k10$membership)
spe$SNN_k10 <- clust_k10 ## Add this one to the SPE too
## Export for later use
cluster_export(
    spe,
    "SNN_k10",
    cluster_dir = file.path(dir_rdata, "clusters_graphbased"),
)

message("Running cut_at() from k = 4 to 28")
clust_k5_list <- lapply(4:28, function(n) {
    message(paste(Sys.time(), "n =", n))
    sort_clusters(igraph::cut_at(g_walk_k10, n = n))
})
names(clust_k5_list) <- paste0("SNN_k10_k", 4:28)

## Add clusters to spe colData
for (i in seq_along(names(clust_k5_list))) {
    colData(spe) <- cbind(colData(spe), clust_k5_list[i])
    ## Add proper name
    colnames(colData(spe))[ncol(colData(spe))] <- names(clust_k5_list)[i]

    ## Export for later use outside the SPE object
    cluster_export(
        spe,
        names(clust_k5_list)[i],
        cluster_dir = file.path(dir_rdata, "clusters_graphbased_cut_at")
    )
}

## make plot
sample_ids <- unique(colData(spe)$sample_id)
pdf(file = file.path(dir_plots, "graph_based_harmony.pdf"))
for (i in seq_along(sample_ids)) {
    for (j in seq_along(names(clust_k5_list))) {
        clus_vals <- unique(clust_k5_list[[j]])
        cols <- Polychrome::palette36.colors(length(clus_vals))
        names(cols) <- clus_vals

        my_plot <- vis_clus(
            spe = spe,
            clustervar = names(clust_k5_list)[j],
            sampleid = sample_ids[i],
            colors = cols,
            auto_crop = FALSE,
            assayname = 'counts',
            ... = paste0(" ", names(clust_k5_list)[j])
        )
        print(my_plot)
    }
}
dev.off()

#   Do offset so we can run BayesSpace. Not here that 'array_row' is not
#   constrained to have max value 77; we instead find the largest 'array_row'
#   value of any sample, and use it to ensure samples are at least 5 rows apart
auto_offset_row <- as.numeric(factor(unique(spe$sample_id))) * (max(spe$array_row) + 5)
names(auto_offset_row) <- unique(spe$sample_id)
spe$row <- spe$array_row + auto_offset_row[spe$sample_id]
spe$col <- spe$array_col

## Save new SPE object
saveRDS(spe, file.path(dir_rdata, "spe_harmony.rds"))

## Object size in GB
## (do this near the end in case lobstr crashes, it's happened to me once)
lobstr::obj_size(spe)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

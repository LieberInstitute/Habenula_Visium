## Required libraries
library("spatialLIBD")
library("here")
library("purrr")
library("harmony")
library("parallel")
library("scater")
library("BiocParallel")
library("ggplot2")
library("ggpubr")
library("scran")
library("Polychrome")
library("sessioninfo")


# Define harmony function that would allow us to specify which reduction to use with SingleCellExperiment object
RunHarmony_mod <- function(
        object,
        group.by.vars,
        #reduction.use = "PCA",
        reduction.use = "GLMPCA_approx", # Top 1000
        dims.use = NULL,
        verbose = TRUE,
        reduction.save = "HARMONY",
        ...) {
    ## Get PCA embeddings
    if (!"PCA" %in% SingleCellExperiment::reducedDimNames(object)) {
        stop("PCA must be computed before running Harmony.")
    }
    ## Get GMLPCA embeddings
    if (!"GLMPCA_approx" %in% SingleCellExperiment::reducedDimNames(object)) {
      stop("GLMPCA must be computed before running Harmony.")
    }
    ## PCA embeddings validations
    pca_embedding <-
        SingleCellExperiment::reducedDim(object, reduction.use)
    if (is.null(dims.use)) {
        dims.use <- seq_len(ncol(pca_embedding))
    }

    if (is.null(dims.use)) {
        dims.use <- seq_len(ncol(pca_embedding))
    }
    dims_avail <- seq_len(ncol(pca_embedding))
    if (!all(dims.use %in% dims_avail)) {
        stop(
            "trying to use more dimensions than computed with PCA. Rerun
            PCA with more dimensions or use fewer PCs"
        )
    }

    metavars_df <- SingleCellExperiment::colData(object)
    if (!all(group.by.vars %in% colnames(metavars_df))) {
        stop("Trying to integrate over variables missing in colData")
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
    colnames(harmonyEmbed) <-
        paste0(reduction.save, "_", seq_len(ncol(harmonyEmbed)))
    SingleCellExperiment::reducedDim(object, reduction.save) <-
        harmonyEmbed

    return(object)
}

tsne_perplex_vals <- c("30", "50") #, "05", "80"
num_cores <- detectCores() - 1

## Create output directories
dir_plots <- here("plots", "04_harmony_BayesSpace")
dir_rdata <- here("processed-data", "04_harmony_BayesSpace")
#harmony_hdf5_dir <- here("processed-data", "04_harmony_BayesSpace", "spe_harmony")

dir.create(dir_plots, showWarnings = FALSE)
dir.create(dir_rdata, showWarnings = FALSE)
dir.create(file.path(dir_rdata, "clusters_graphbased"), showWarnings = FALSE)
dir.create(file.path(dir_rdata, "clusters_graphbased_cut_at"), showWarnings = FALSE)

set.seed(20240614)

## Load the data
# spe <- loadHDF5SummarizedExperiment(filtered_hdf5_dir)
spe <- readRDS(file.path(dir_rdata, "spe_qcED_spatialLIBD_log_GLM-PCA.rds")) 
# colnames(colData(spe))
reducedDimNames(spe)
# [1] "10x_pca"       "10x_tsne"      "10x_umap"      "PCA"          
# [5] "GLMPCA_approx"

## Plot initial low-dimensional representations prior to batch correction

# Build a list with the reductions PCA and GLM-PCA to plot

lst_PCA <- c(reducedDimNames(spe)[grep("^PCA", reducedDimNames(spe))])

## DimRed plots before Harmony
pdf(file.path(dir_plots, 'reduction_dimension_PCA.pdf'), useDingbats = FALSE)
map(lst_PCA, ~ plotReducedDim(spe,
                              dimred = .x, 
                              ncomponents = 3,
                              colour_by = "sample_id"
)) 
dev.off()

lst_GLMPCA <- c(reducedDimNames(spe)[grep("^GLMPCA_", reducedDimNames(spe))])
pdf(file.path(dir_plots, 'reduction_dimension_GLMPCA.pdf'), useDingbats = FALSE)
map(lst_GLMPCA, ~ plotReducedDim(spe,
                              dimred = .x, 
                              ncomponents = 3,
                              colour_by = "sample_id"
)) 
dev.off()

## plot other GML-PCA features
point_size = 0.5
font_size = 7
pdf(file.path(dir_plots, 'reduction_dimension_GLMPCA_other_features.pdf'), useDingbats = FALSE)
plotReducedDim(spe, dimred = "GLMPCA_approx", ncomponents = 3, colour_by = "brain_id", point_size = point_size) + 
  theme(axis.text=element_text(size=font_size) ,axis.title=element_text(size=font_size))
plotReducedDim(spe, dimred = "GLMPCA_approx", ncomponents = 3, colour_by = "brain_area", point_size = point_size) + 
  theme(axis.text=element_text(size=font_size) ,axis.title=element_text(size=font_size))
plotReducedDim(spe, dimred = "GLMPCA_approx", ncomponents = 3, colour_by = "ethnicity", point_size = point_size) + 
  theme(axis.text=element_text(size=font_size) ,axis.title=element_text(size=font_size))
plotReducedDim(spe, dimred = "GLMPCA_approx", ncomponents = 3, colour_by = "sex", point_size = point_size) + 
  theme(axis.text=element_text(size=font_size) ,axis.title=element_text(size=font_size))
#gridExtra::grid.arrange(plt1, plt2, plt3, plt4, ncol=1, nrow=2)
dev.off()

## Perform harmony batch correction
message("Running RunHarmony()")
Sys.time()
set.seed(20240614)

## add additional co-variables
colnames(colData(spe))
#covars <- c("sample_id") # This is the highest technical level / includes subsets of "brain_id" and "ethnicity"
covars <- c("sample_id", "brain_id")

spe <-
    RunHarmony_mod(
        spe,
        group.by.vars = covars, #"sample_id",
        verbose = TRUE,
        plot_convergence = TRUE,
        #reduction.use = "PCA",   # HVGs at FDR = 0.05 = 8009
        reduction.use = "GLMPCA_approx",  # Top 1000
        reduction.save = "HARMONY",
        kmeans_init_nstart = 100,
        kmeans_init_iter_max = 1000
    )

# spe <-
#     RunHarmony_mod(
#         spe,
#         group.by.vars = covars, #"sample_id",
#         verbose = TRUE,
#         #reduction.use = "PCA",
#         reduction.use = "GLMPCA_approx",  # Top 1000
#         reduction.save = "harmony_subject_no_lambda",
#         plot_convergence = TRUE,
#         lambda = NULL,
#         max_iter = 30
#     )

Sys.time()
reducedDimNames(spe)
# [1] "10x_pca"       "10x_tsne"      "10x_umap"      "PCA"          
# [5] "GLMPCA_approx" "HARMONY"     

## Run Harmony on GLMPCA too, with and without lambda = NULL

#   Perform dimensionality reduction using both PCA and harmony's reduced
#   dimensions
# for (dimred_var in c("PCA", "HARMONY", "harmony_subject_no_lambda")) {
for (dimred_var in c("GLMPCA_approx", "HARMONY")) {
    #   Run TSNE with several perplexity values
    for (perplex in tsne_perplex_vals) {
        message(
            sprintf(
                "Running runTSNE() perplexity %s on %s dimensions",
                perplex,
                dimred_var
            )
        )
        Sys.time()
        set.seed(20240229)
        spe <-
            runTSNE(
                spe,
                dimred = dimred_var,
                name = sprintf("TSNE_perplexity%s.%s", perplex, dimred_var),
                perplexity = as.integer(perplex)
            )
        Sys.time()

        #   Explore TSNE results via plots
        p_tsne <- ggplot(
            data.frame(reducedDim(
                spe,
                sprintf("TSNE_perplexity%s.%s", perplex, dimred_var)
            )),
            aes(
                x = TSNE1,
                y = TSNE2,
                color = factor(spe$sample_id)
            )
        ) +
            geom_point() +
            labs(color = "sample_id") +
            theme_bw()
        pdf(
            file = file.path(
                dir_plots,
                sprintf(
                    "tSNE_perplexity%s_%s_sample_id.pdf",
                    perplex,
                    dimred_var
                )
            ),
            width = 9
        )
        print(p_tsne)
        dev.off()
    }

    #   Also run UMAP
    message(sprintf("Running runUMAP() on %s dimensions", dimred_var))
    Sys.time()
    set.seed(20240229)
    spe <- runUMAP(
        spe,
        dimred = dimred_var,
        name = sprintf("UMAP.%s", dimred_var),
        BPPARAM = MulticoreParam(num_cores)
    )
    Sys.time()

    #   Explore UMAP results, coloring by both subject and sample ID
    for (color_var in c("sample_id")) {
        p_umap <- ggplot(
            data.frame(reducedDim(
                spe, sprintf("UMAP.%s", dimred_var)
            )),
            aes(
                x = UMAP1,
                y = UMAP2,
                color = factor(spe[[color_var]])
            )
        ) +
            geom_point() +
            labs(color = color_var) +
            theme_bw()

        pdf(file = file.path(
            dir_plots,
            sprintf("UMAP_%s_%s.pdf", color_var, dimred_var)
        ))
        print(p_umap)
        dev.off()
    }
}

##  Plot the UMAP by sum of UMIs to verified low/zero UMIs
colnames(colData(spe))
reducedDimNames(spe)

## plots after Harmony to compare against GLMPCA_approx

plotsGLMPCA_Harmony_facets <- function(dimred_name, feature_by) {

  point_size = 0.5
  point_sizef = 0.3
  title_name <- paste0(dimred_name, ": ", feature_by)
  
  plt1 <- plotReducedDim(spe, dimred = dimred_name, colour_by = feature_by, point_size = point_size) + ggtitle(title_name) 
  plt1b <- plotReducedDim(spe, dimred = dimred_name, colour_by = feature_by, point_size = point_sizef) +
    facet_wrap(~ spe[[feature_by]]) +
    theme(legend.position="none",
          plot.background = element_rect(fill = "white", color = "white"),  # White background
          panel.background = element_rect(fill = "white", color = "white"))  # White panel bg
  plt1_f <- ggarrange(plt1, plt1b, ncol = 1, nrow = 2)

  return(plt1_f)
}


  
pdf(file.path(dir_plots, 'reduction_dimension_Harmony_vs_GLMPCA.pdf'), useDingbats = FALSE)

plt1 <- plotsGLMPCA_Harmony_facets("GLMPCA_approx", "sample_id")
plt2 <- plotsGLMPCA_Harmony_facets("HARMONY", "sample_id")
plt3 <- plotsGLMPCA_Harmony_facets("GLMPCA_approx", "brain_id")
plt4 <- plotsGLMPCA_Harmony_facets("HARMONY", "brain_id")
plt5 <- plotsGLMPCA_Harmony_facets("GLMPCA_approx", "brain_area")
plt6 <- plotsGLMPCA_Harmony_facets("HARMONY", "brain_area") 
plt7 <- plotsGLMPCA_Harmony_facets("GLMPCA_approx", "ethnicity")
plt8 <- plotsGLMPCA_Harmony_facets("HARMONY", "ethnicity")

# plt1 <- plotReducedDim(spe, dimred = "GLMPCA_approx", colour_by = "sample_id", point_size = point_size)
# plt2 <- plotReducedDim(spe, dimred = "HARMONY", colour_by = "sample_id", point_size = point_size)
# plt3 <- plotReducedDim(spe, dimred = "GLMPCA_approx", colour_by = "brain_id", point_size = point_size) 
# plt4 <- plotReducedDim(spe, dimred = "HARMONY", colour_by = "brain_id", point_size = point_size)
# plt5 <- plotsGLMPCA_Harmony_facets("GLMPCA_approx", "brain_area", "GLMPCA: brain_area")
# plt6 <- plotReducedDim(spe, dimred = "HARMONY", colour_by = "brain_area", point_size = point_size) 
# plt7 <- plotReducedDim(spe, dimred = "GLMPCA_approx", colour_by = "ethnicity", point_size = point_size) 
# plt8 <- plotReducedDim(spe, dimred = "HARMONY", colour_by = "ethnicity", point_size = point_size)

plt_all <- ggarrange(plt1, plt2, plt3, plt4, plt5, plt6, plt7, plt8, 
                     labels = c("A", "B", "C", "D"),
                     ncol = 1, nrow = 1)
plt_all
dev.off()

## plots sum_umi from Harmony against GLMPCA_approx colored by brain_id (3 levels)
# Note shape_by is restricted to 10 levels
pdf(file.path(dir_plots, 'reduction_dimension_Harmony_vs_GLMPCA_sumUMI_sumGene.pdf'), useDingbats = FALSE)
plt1 <- plotReducedDim(spe, dimred = "HARMONY", by_exprs_values = "logcounts", shape_by = "brain_id", colour_by = "sum_umi") + 
  ggtitle("HARMONY (GLMPCA - 1000 HDVG)")
plt2 <- plotReducedDim(spe, dimred = "HARMONY", by_exprs_values = "logcounts", shape_by = "brain_id", colour_by = "sum_gene")
plt3 <- plotReducedDim(spe, dimred = "GLMPCA_approx", by_exprs_values = "logcounts", shape_by = "brain_id", colour_by = "sum_umi") + 
  ggtitle("GLMPCA (1000 HDVG)") 
plt4 <- plotReducedDim(spe, dimred = "GLMPCA_approx", by_exprs_values = "logcounts", shape_by = "brain_id", colour_by = "sum_gene")
plt_all <- ggarrange(plt1, plt2, plt3, plt4,
                        labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2, common.legend = TRUE, legend="right")
plt_all
dev.off()


## Perform graph-based clustering on batch corrected-data. Smaller 'k' usually yields finer clusters (ex. 5)
message("Running buildSNNGraph() on HARMONY dimensions")
Sys.time()
g_k10 <- buildSNNGraph(spe, k = 10, use.dimred = "HARMONY")
Sys.time()
save(g_k10, file = file.path(dir_rdata, "g_k10_harmony.Rdata"))

## For clustering based on the produced graph
message("Running cluster_walktrap()")
Sys.time()
g_walk_k10 <- igraph::cluster_walktrap(g_k10)
Sys.time()
save(g_walk_k10, file = file.path(dir_rdata, "g_walk_k10_harmony.Rdata"))

clust_k10 <- sort_clusters(g_walk_k10$membership)
spe$SNN_k10 <- clust_k10 ## Add this one to the SPE too
## Export for later use
cluster_export(spe,
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
    colnames(colData(spe))[ncol(colData(spe))] <-
        names(clust_k5_list)[i]

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
            assayname = "counts",
            ... = paste0(" ", names(clust_k5_list)[j])
        )
        print(my_plot)
    }
}
dev.off()

#   Do offset so we can run BayesSpace. Not here that 'array_row' is not
#   constrained to have max value 77; we instead find the largest 'array_row'
#   value of any sample, and use it to ensure samples are at least 5 rows apart
auto_offset_row <-
    as.numeric(factor(unique(spe$sample_id))) * (max(spe$array_row) + 5)
names(auto_offset_row) <- unique(spe$sample_id)
spe$row <- spe$array_row + auto_offset_row[spe$sample_id]
spe$col <- spe$array_col

## Save new SPE object
saveRDS(spe, file.path(dir_rdata, "spe_harmony.rds"))

## Object size in GB. 4.44 GB
## (do this near the end in case lobstr crashes, it's happened to me once)
lobstr::obj_size(spe)

message('Harmony correction completed! ')   


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# > print("Reproducibility information:")
# [1] "Reproducibility information:"
# > Sys.time()
# [1] "2024-06-14 12:46:38 EDT"
# > proc.time()
# user   system  elapsed 
# 1536.347   14.137 7150.028 
# > options(width = 120)
# > session_info()
# rkmeData          1.0.4       2020-04-23 [2] CRAN (R 4.3.2)
# Biobase                * 2.62.0      2023-10-24 [2] Bioconductor
# BiocFileCache            2.10.1      2023-10-26 [2] Bioconductor
# BiocGenerics           * 0.48.1      2023-11-01 [2] Bioconductor
# BiocIO                   1.12.0      2023-10-24 [2] Bioconductor
# BiocManager              1.30.22     2023-08-08 [2] CRAN (R 4.3.2)
# BiocNeighbors            1.20.2      2024-01-07 [2] Bioconductor 3.18 (R 4.3.2)
# BiocParallel           * 1.36.0      2023-10-24 [2] Bioconductor
# BiocSingular             1.18.0      2023-10-24 [2] Bioconductor
# BiocVersion              3.18.1      2023-11-15 [2] Bioconductor
# Biostrings               2.70.2      2024-01-28 [2] Bioconductor 3.18 (R 4.3.2)
# bit                      4.0.5       2022-11-15 [2] CRAN (R 4.3.2)
# bit64                    4.0.5       2020-08-30 [2] CRAN (R 4.3.2)
# bitops                   1.0-7       2021-04-24 [2] CRAN (R 4.3.2)
# blob                     1.2.4       2023-03-17 [2] CRAN (R 4.3.2)
# bluster                  1.12.0      2023-10-24 [2] Bioconductor
# bslib                    0.6.1       2023-11-28 [2] CRAN (R 4.3.2)
# cachem                   1.0.8       2023-05-01 [2] CRAN (R 4.3.2)
# cli                      3.6.2       2023-12-11 [2] CRAN (R 4.3.2)
# cluster                  2.1.6       2023-12-01 [3] CRAN (R 4.3.2)
# codetools                0.2-19      2023-02-01 [3] CRAN (R 4.3.2)
# colorspace               2.1-0       2023-01-23 [2] CRAN (R 4.3.2)
# config                   0.3.2       2023-08-30 [2] CRAN (R 4.3.2)
# cowplot                  1.1.3       2024-01-22 [2] CRAN (R 4.3.2)
# crayon                   1.5.2       2022-09-29 [2] CRAN (R 4.3.2)
# curl                     5.2.0       2023-12-08 [2] CRAN (R 4.3.2)
# data.table               1.15.0      2024-01-30 [2] CRAN (R 4.3.2)
# DBI                      1.2.1       2024-01-12 [2] CRAN (R 4.3.2)
# dbplyr                   2.4.0       2023-10-26 [2] CRAN (R 4.3.2)
# DelayedArray             0.28.0      2023-10-24 [2] Bioconductor
# DelayedMatrixStats       1.24.0      2023-10-24 [2] Bioconductor
# digest                   0.6.34      2024-01-11 [2] CRAN (R 4.3.2)
# doParallel               1.0.17      2022-02-07 [2] CRAN (R 4.3.2)
# dotCall64                1.1-1       2023-11-28 [2] CRAN (R 4.3.2)
# dplyr                    1.1.4       2023-11-17 [2] CRAN (R 4.3.2)
# dqrng                    0.3.2       2023-11-29 [2] CRAN (R 4.3.2)
# DT                       0.31        2023-12-09 [2] CRAN (R 4.3.2)
# edgeR                    4.0.14      2024-01-29 [2] Bioconductor 3.18 (R 4.3.2)
# ellipsis                 0.3.2       2021-04-29 [2] CRAN (R 4.3.2)
# ExperimentHub            2.10.0      2023-10-24 [2] Bioconductor
# fansi                    1.0.6       2023-12-08 [2] CRAN (R 4.3.2)
# farver                   2.1.1       2022-07-06 [2] CRAN (R 4.3.2)
# fastmap                  1.1.1       2023-02-24 [2] CRAN (R 4.3.2)
# fields                   15.2        2023-08-17 [2] CRAN (R 4.3.2)
# filelock                 1.0.3       2023-12-11 [2] CRAN (R 4.3.2)
# foreach                  1.5.2       2022-02-02 [2] CRAN (R 4.3.2)
# generics                 0.1.3       2022-07-05 [2] CRAN (R 4.3.2)
# GenomeInfoDb           * 1.38.5      2023-12-28 [2] Bioconductor 3.18 (R 4.3.2)
# GenomeInfoDbData         1.2.11      2024-02-09 [2] Bioconductor
# GenomicAlignments        1.38.2      2024-01-16 [2] Bioconductor 3.18 (R 4.3.2)
# GenomicRanges          * 1.54.1      2023-10-29 [2] Bioconductor
# ggbeeswarm               0.7.2       2023-04-29 [2] CRAN (R 4.3.2)
# ggplot2                * 3.4.4       2023-10-12 [2] CRAN (R 4.3.2)
# ggrepel                  0.9.5       2024-01-10 [2] CRAN (R 4.3.2)
# glue                     1.7.0       2024-01-09 [2] CRAN (R 4.3.2)
# golem                    0.4.1       2023-06-05 [2] CRAN (R 4.3.2)
# gridExtra                2.3         2017-09-09 [2] CRAN (R 4.3.2)
# gtable                   0.3.4       2023-08-21 [2] CRAN (R 4.3.2)
# harmony                * 1.2.0       2023-11-29 [2] CRAN (R 4.3.2)
# here                   * 1.0.1       2020-12-13 [2] CRAN (R 4.3.2)
# htmltools                0.5.7       2023-11-03 [2] CRAN (R 4.3.2)
# htmlwidgets              1.6.4       2023-12-06 [2] CRAN (R 4.3.2)
# httpuv                   1.6.14      2024-01-26 [2] CRAN (R 4.3.2)
# httr                     1.4.7       2023-08-15 [2] CRAN (R 4.3.2)
# igraph                   2.0.1.9008  2024-02-09 [2] Github (igraph/rigraph@39158c6)
# interactiveDisplayBase   1.40.0      2023-10-24 [2] Bioconductor
# IRanges                * 2.36.0      2023-10-24 [2] Bioconductor
# irlba                    2.3.5.1     2022-10-03 [2] CRAN (R 4.3.2)
# iterators                1.0.14      2022-02-05 [2] CRAN (R 4.3.2)
# jquerylib                0.1.4       2021-04-26 [2] CRAN (R 4.3.2)
# jsonlite                 1.8.8       2023-12-04 [2] CRAN (R 4.3.2)
# KEGGREST                 1.42.0      2023-10-24 [2] Bioconductor
# labeling                 0.4.3       2023-08-29 [2] CRAN (R 4.3.2)
# later                    1.3.2       2023-12-06 [2] CRAN (R 4.3.2)
# lattice                  0.22-5      2023-10-24 [3] CRAN (R 4.3.2)
# lazyeval                 0.2.2       2019-03-15 [2] CRAN (R 4.3.2)
# lifecycle                1.0.4       2023-11-07 [2] CRAN (R 4.3.2)
# limma                    3.58.1      2023-10-31 [2] Bioconductor
# locfit                   1.5-9.8     2023-06-11 [2] CRAN (R 4.3.2)
# magick                   2.8.2       2023-12-20 [2] CRAN (R 4.3.2)
# magrittr                 2.0.3       2022-03-30 [2] CRAN (R 4.3.2)
# maps                     3.4.2       2023-12-15 [2] CRAN (R 4.3.2)
# Matrix                   1.6-5       2024-01-11 [3] CRAN (R 4.3.2)
# MatrixGenerics         * 1.14.0      2023-10-24 [2] Bioconductor
# matrixStats            * 1.2.0       2023-12-11 [2] CRAN (R 4.3.2)
# memoise                  2.0.1       2021-11-26 [2] CRAN (R 4.3.2)
# metapod                  1.10.1      2023-12-24 [2] Bioconductor 3.18 (R 4.3.2)
# mime                     0.12        2021-09-28 [2] CRAN (R 4.3.2)
# munsell                  0.5.0       2018-06-12 [2] CRAN (R 4.3.2)
# paletteer                1.6.0       2024-01-21 [2] CRAN (R 4.3.2)
# pillar                   1.9.0       2023-03-22 [2] CRAN (R 4.3.2)
# pkgconfig                2.0.3       2019-09-22 [2] CRAN (R 4.3.2)
# plotly                   4.10.4      2024-01-13 [2] CRAN (R 4.3.2)
# png                      0.1-8       2022-11-29 [2] CRAN (R 4.3.2)
# Polychrome             * 1.5.1       2022-04-29 [1] R-Forge (R 4.3.2)
# promises                 1.2.1       2023-08-10 [2] CRAN (R 4.3.2)
# purrr                  * 1.0.2       2023-08-10 [2] CRAN (R 4.3.2)
# R6                       2.5.1       2021-08-19 [2] CRAN (R 4.3.2)
# rappdirs                 0.3.3       2021-01-31 [2] CRAN (R 4.3.2)
# RColorBrewer             1.1-3       2022-04-03 [2] CRAN (R 4.3.2)
# Rcpp                   * 1.0.12      2024-01-09 [2] CRAN (R 4.3.2)
# RcppAnnoy                0.0.22      2024-01-23 [2] CRAN (R 4.3.2)
# RCurl                    1.98-1.14   2024-01-09 [2] CRAN (R 4.3.2)
# rematch2                 2.1.2       2020-05-01 [2] CRAN (R 4.3.2)
# restfulr                 0.0.15      2022-06-16 [2] CRAN (R 4.3.2)
# RhpcBLASctl              0.23-42     2023-02-11 [2] CRAN (R 4.3.2)
# rjson                    0.2.21      2022-01-09 [2] CRAN (R 4.3.2)
# rlang                    1.1.3       2024-01-10 [2] CRAN (R 4.3.2)
# rprojroot                2.0.4       2023-11-05 [2] CRAN (R 4.3.2)
# Rsamtools                2.18.0      2023-10-24 [2] Bioconductor
# RSQLite                  2.3.5       2024-01-21 [2] CRAN (R 4.3.2)
# rsvd                     1.0.5       2021-04-16 [2] CRAN (R 4.3.2)
# rtracklayer              1.62.0      2023-10-24 [2] Bioconductor
# Rtsne                    0.17        2023-12-07 [2] CRAN (R 4.3.2)
# S4Arrays                 1.2.0       2023-10-24 [2] Bioconductor
# S4Vectors              * 0.40.2      2023-11-23 [2] Bioconductor 3.18 (R 4.3.2)
# sass                     0.4.8       2023-12-06 [2] CRAN (R 4.3.2)
# ScaledMatrix             1.10.0      2023-10-24 [2] Bioconductor
# scales                   1.3.0       2023-11-28 [2] CRAN (R 4.3.2)
# scater                 * 1.30.1      2023-11-16 [2] Bioconductor
# scatterplot3d            0.3-44      2023-05-05 [1] R-Forge (R 4.3.2)
# scran                  * 1.30.2      2024-01-22 [2] Bioconductor 3.18 (R 4.3.2)
# scuttle                * 1.12.0      2023-10-24 [2] Bioconductor
# sessioninfo            * 1.2.2       2021-12-06 [2] CRAN (R 4.3.2)
# shiny                    1.8.0       2023-11-17 [2] CRAN (R 4.3.2)
# shinyWidgets             0.8.1       2024-01-10 [2] CRAN (R 4.3.2)
# SingleCellExperiment   * 1.24.0      2023-10-24 [2] Bioconductor
# spam                     2.10-0      2023-10-23 [2] CRAN (R 4.3.2)
# SparseArray              1.2.3       2023-12-25 [2] Bioconductor 3.18 (R 4.3.2)
# sparseMatrixStats        1.14.0      2023-10-24 [2] Bioconductor
# SpatialExperiment      * 1.12.0      2023-10-24 [2] Bioconductor
# spatialLIBD            * 1.15.4      2024-05-01 [1] Github (LieberInstitute/spatialLIBD@77a5303)
# statmod                  1.5.0       2023-01-06 [2] CRAN (R 4.3.2)
# SummarizedExperiment   * 1.32.0      2023-10-24 [2] Bioconductor
# tibble                   3.2.1       2023-03-20 [2] CRAN (R 4.3.2)
# tidyr                    1.3.1       2024-01-24 [2] CRAN (R 4.3.2)
# tidyselect               1.2.0       2022-10-10 [2] CRAN (R 4.3.2)
# utf8                     1.2.4       2023-10-22 [2] CRAN (R 4.3.2)
# uwot                     0.1.16      2023-06-29 [2] CRAN (R 4.3.2)
# vctrs                    0.6.5       2023-12-01 [2] CRAN (R 4.3.2)
# vipor                    0.4.7       2023-12-18 [2] CRAN (R 4.3.2)
# viridis                  0.6.5       2024-01-29 [2] CRAN (R 4.3.2)
# viridisLite              0.4.2       2023-05-02 [2] CRAN (R 4.3.2)
# withr                    3.0.0       2024-01-16 [2] CRAN (R 4.3.2)
# XML                      3.99-0.16.1 2024-01-22 [2] CRAN (R 4.3.2)
# xtable                   1.8-4       2019-04-21 [2] CRAN (R 4.3.2)
# XVector                  0.42.0      2023-10-24 [2] Bioconductor
# yaml                     2.3.8       2023-12-11 [2] CRAN (R 4.3.2)
# zlibbioc                 1.48.0      2023-10-24 [2] Bioconductor
# 
# [1] /users/csoto/R/4.3.x
# [2] /jhpce/shared/community/core/conda_R/4.3.x/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.3.x/R/lib64/R/library

library("spatialLIBD")
library("scran")
library("tidyverse")
library("dplyr")
library("here")
library("lobstr")
library("sessioninfo")


## Create directory plots

dir_plots <- here::here("plots", "02_build_spe")
if (!dir.exists(dir_plots)) {
    dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
}

## Set directory data

dir_rdata <- here("processed-data", "02_build_spe")

## set path to raw and pre-filtered data

spe_in_path <- here("processed-data", "02_build_spe", "spe.rds")
raw_in_path <- here("processed-data", "02_build_spe", "spe_raw.rds")

# /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/02_build_spe/spe.rds


## load Datasets

spe_raw <- readRDS(raw_in_path)
message("Initial number of spots:", dim(spe_raw)[2], "\n")
# merged samples: Initial number of spots: 44928

spe <- readRDS(spe_in_path)
message("Initial number of spots :", dim(spe)[2], "\n")
# merged samples: Initial number of spots: 16928

## Set some initials for manage spot size in the plots
var_height <- 24 # 24/3=8
var_width <- 26 # 36/4=9
var_point_size <- 1.5

set.seed(07112024)

# origin code copied from https://github.com/LieberInstitute/spatialDLPFC/blob/14a1f253a92e43c01fec3cc3077a9b2cf9ce9fc0/code/analysis/01_build_spe/01_build_spe.R#L216-L221

## Inspect in vs outside of tissue

## Get number of TRUE spots in tissue
in_tissue_spots <- map(unique(spe_raw$sample_id), ~ summary(spe_raw$in_tissue[spe_raw$sample_id == .x] == TRUE))
in_tissue_spots_F <- sum(as.numeric(sapply(in_tissue_spots, "[[", 2)))
in_tissue_spots_T <- sum(as.numeric(sapply(in_tissue_spots, "[[", 3)))

print(paste0("Spots in tissue FALSE: ", in_tissue_spots_F, " TRUE: ", in_tissue_spots_T))

lst_order <- sort(unique(spe$sample_id))

vis_grid_clus(
    spe = spe_raw,
    clustervar = "in_tissue",
    sample_order = lst_order,
    height = var_height, # 8
    width = var_width, # 9
    point_size = var_point_size,
    pdf = here(dir_plots, "all_in_tissue_grid.pdf"),
    sort_clust = FALSE,
    colors = c("TRUE" = "grey90", "FALSE" = "orange")
)


## -----------------------------
## Out-tissue metrics

lst_out_counts <- c(
    sum_umi = "out_tissue_sum_umi.pdf",
    sum_gene = "out_tissue_sum_gene.pdf",
    expr_chrM_ratio = "out_tissue_expr_chrM_ratio.pdf"
)

print("Ploting out-tissues metrics")

map2(as.vector(names(lst_out_counts)), as.vector(lst_out_counts), ~ vis_grid_gene(
    spe = spe_raw[, which(!colData(spe_raw)$in_tissue)],
    geneid = .x,
    sample_order = lst_order,
    height = var_height, # 8
    width = var_width, # 9
    point_size = var_point_size,
    # return_plots = TRUE,
    pdf = here(dir_plots, .y),
    assayname = "counts"
))


print("Plots done!")

summary(spe_raw$sum_umi[which(!colData(spe_raw)$in_tissue)])
mean(spe_raw$sum_umi[which(!colData(spe_raw)$in_tissue)])

unique(spe_raw$sample_id)
map(unique(spe_raw$sample_id), ~ summary(spe_raw$sum_umi[spe_raw$sample_id == .x]))
map(unique(spe_raw$sample_id), ~ median(spe_raw$sum_umi[spe_raw$sample_id == .x]))

summary(spe_raw$sum_gene[which(!colData(spe_raw)$in_tissue)])

map(unique(spe_raw$sample_id), ~ summary(spe_raw$sum_gene[spe_raw$sample_id == .x]))

summary(spe_raw$expr_chrM_ratio[which(!colData(spe_raw)$in_tissue)])

map(unique(spe_raw$sample_id), ~ summary(spe_raw$expr_chrM_ratio[spe_raw$sample_id == .x]))

# head(table(spe_raw$sum_umi[which(!colData(spe_raw)$in_tissue)]))
# # 202 428 432 436 455 482
# # 1   1   1   1   1   1

# Metrics measured on the SpaceRanger report correspond to the spe object; e.g. This is the metric called "Median UMI counts per spot"
map(unique(spe$sample_id), ~ median(spe$sum_umi[spe$sample_id == .x])) 
map(unique(spe$sample_id), ~ summary(spe$sum_gene[spe$sample_id == .x]))

## -----------------------------
## in-tissue metrics

lst_in_counts <- c(
    sum_umi = "in_tissue_sum_umi.pdf",
    sum_gene = "in_tissue_sum_gene.pdf",
    expr_chrM_ratio = "in_tissue_expr_chrM_ratio.pdf"
)

print("Ploting in-tissues plots")

map2(as.vector(names(lst_in_counts)), as.vector(lst_in_counts), ~ vis_grid_gene(
    spe = spe,
    geneid = .x,
    sample_order = lst_order,
    height = var_height, # 8
    width = var_width, # 9
    point_size = var_point_size,
    # return_plots = TRUE,
    pdf = here(dir_plots, .y),
    assayname = "counts")
    )


## Calculate total genes with count the first few gene sums
counts_matrix <- (assay(spe, "counts"))
## Combine gene names and their sums
sum_genes <- rowSums(assay(spe, "counts"))
## Add these sums to the rowData of the spe object
rowData(spe)$sum_counts <- sum_genes
gene_summary <- data.frame(
  gene_name_ens = rownames(spe),
  gene_name = rowData(spe)["gene_name"],
  sum_counts = sum_genes
)
gene_summary <- gene_summary[order(gene_summary$sum_counts, decreasing = F), ]
tail(gene_summary)
#                    gene_name_ens gene_name sum_counts
# ENSG00000198804 ENSG00000198804    MT-CO1    3656352
# ENSG00000198712 ENSG00000198712    MT-CO2    2934506
f_name <- paste0(here(dir_rdata, "spe_gene_counts.csv"))
write.csv(gene_summary, f_name, row.names = FALSE)

# Group UMI sums by slide
# Sum UMIs per spot (column-wise sum)
spot_UMI_sums <- colSums(counts_matrix)
# Add UMI sums to colData
colData(spe)$UMI_sum <- spot_UMI_sums
colnames(colData(spe))
UMI_per_slide <- colData(spe) %>%
  as.data.frame() %>%
  group_by(sample_id) %>%
  summarise(total_UMIs = sum(UMI_sum))
head(UMI_per_slide)

# Count spots per slide
spots_by_slide <- table(colData(spe)$sample_id)
# Convert to a data frame for easier manipulation if needed
spots_by_slide_df <- as.data.frame(spots_by_slide)
colnames(spots_by_slide_df) <- c("Slide", "Total_Spots")
#print(spots_by_slide_df)
cbind(spots_by_slide_df, UMI_per_slide)

summary(spe$sum_umi)

map(unique(spe$sample_id), ~ summary(spe$sum_umi[spe$sample_id == .x]))

summary(spe$sum_gene)

map(unique(spe$sample_id), ~ summary(spe$sum_gene[spe$sample_id == .x]))

summary(spe$expr_chrM_ratio)

map(unique(spe$sample_id), ~ summary(spe$expr_chrM_ratio[spe$sample_id == .x]))



## -----------------------------
## All in and out tissue metrics

lst_all_counts <- c(
    sum_umi = "all_sum_umi.pdf",
    sum_gene = "all_sum_gene.pdf",
    expr_chrM_ratio = "all_expr_chrM_ratio.pdf"
)

print("Ploting ALL in and out tissues plots")

map2(as.vector(names(lst_all_counts)), as.vector(lst_all_counts), ~ vis_grid_gene(
    spe = spe_raw,
    geneid = .x,
    sample_order = lst_order,
    height = var_height, # 8
    width = var_width, # 9
    point_size = var_point_size,
    # return_plots = TRUE,
    pdf = here(dir_plots, .y),
    assayname = "counts"
))


## -----------------------------
## Add QC metrics, followed by filtering data on edge spots

## Adapting code from https://github.com/LieberInstitute/Habenula_Visium/blob/0e020dd1a580b2bcea130057a2d7c01e5f928001/code/04_harmony_BayesSpace/01-filter_normalize.R#L29C1-L150C1

message("Initial number of spots:", dim(spe)[2], "\n")
# Preliminary QC
spe <- spe[
    rowSums(assays(spe)$counts) > 0,
    (colSums(assays(spe)$counts) > 0) & spe$in_tissue]

message("Number of spots after preliminary QC:", dim(spe)[2], "\n")
# Number of spots after preliminary QC: 16928


## Metrics QC

metrics_qc <- function(spe) {
  
  spe$in_tissue <- as.logical(spe$in_tissue)
  spe_in <- spe[, spe$in_tissue]
  
  ## QC in-tissue spots
  ## define variables  
  
  qc_df <- data.frame(
        log2sum = log2(spe$sum_umi),
        log2detected = log2(spe$sum_gene),
        subsets_Mito_percent = spe$expr_chrM_ratio * 100,
        sample_id = spe$sample_id)
    #head(qc_df)
  
   qcfilter <- data.frame(
        low_lib_size = scater::isOutlier(qc_df$log2sum, type = "lower", log = FALSE, batch = qc_df$sample_id),
        low_n_features = scater::isOutlier(qc_df$log2detected, type = "lower", log = FALSE, batch = qc_df$sample_id),
        high_subsets_Mito_percent = scater::isOutlier(qc_df$subsets_Mito_percent, type = "higher", batch = qc_df$sample_id)
    ) |>
      dplyr::mutate(discard = (low_lib_size | low_n_features) | high_subsets_Mito_percent)

    ## Add qc filter cols to colData(spe) after factoring
    ## low_lib_size_low_mito / discard / low_lib_size / low_n_features / high mito percent
    
    spe$scran_low_lib_size_low_mito <- factor(qcfilter$low_lib_size & qc_df$subsets_Mito_percent < 0.5, levels = c("TRUE", "FALSE"))
    spe$scran_discard <- factor(qcfilter$discard, levels = c("TRUE", "FALSE"))
    spe$scran_low_lib_size <- factor(qcfilter$low_lib_size, levels = c("TRUE", "FALSE"))
    spe$scran_low_n_features <- factor(qcfilter$low_n_features, levels = c("TRUE", "FALSE"))
    spe$scran_high_subsets_Mito_percent <- factor(qcfilter$high_subsets_Mito_percent, levels = c("TRUE", "FALSE"))

    ## Find edge spots
    
    array_row <- array_col <- edge_row <- edge_col <- row_distance <- NULL
    col_distance <- high_subsets_Mito_percent <- NULL
    
    spot_coords <- colData(spe_in) |>
      as.data.frame() |>
      select(in_tissue, sample_id, array_row, array_col) |>
      group_by(sample_id, array_row) |>
      mutate(
        edge_col = array_col == min(array_col) | array_col == max(array_col),
        col_distance = pmin(
          abs(array_col - min(array_col)),
          abs(array_col - max(array_col))
        )
      ) |>
      group_by(sample_id, array_col) |>
      mutate(
        edge_row = array_row == min(array_row) | array_row == max(array_row),
        row_distance = pmin(
          abs(array_row - min(array_row)),
          abs(array_row - max(array_row))
        )
      ) |>
      group_by(sample_id) |>
      mutate(
        edge_spot = edge_row | edge_col,
        edge_distance = pmin(row_distance, col_distance)
      )
    
    
    ## Add Edge info to spe
    spe$edge_spot <- NA
    spe$edge_spot[which(spe$in_tissue)] <- spot_coords$edge_spot
    
    spe$edge_distance <- NA
    spe$edge_distance[which(spe$in_tissue)] <- spot_coords$edge_distance
    
    spe$scran_low_lib_size_edge <- NA
    spe$scran_low_lib_size_edge[which(spe$in_tissue)] <- qcfilter$low_lib_size & spot_coords$edge_spot
    
    # spots <- data.frame(
    #     row = spe$array_row,
    #     col = spe$array_col,
    #     sample_id = spe$sample_id
    # )
    # 
    # edge_spots_row <- group_by(spots, sample_id, row) %>% mutate(min_col = min(col), max_col = max(col))
    # edge_spots_col <- group_by(spots, sample_id, col) %>% mutate(min_row = min(row), max_row = max(row))
    # spots <- left_join(spots, edge_spots_row) %>% left_join(edge_spots_col)
    # spots$edge_spots <-
    #     with(
    #         spots,
    #         row == min_row | row == max_row | col == min_col | col == max_col
    #     )
    # 
    # head(spots, n = 3)
    # # row col     sample_id min_col max_col min_row max_row edge_spots
    # # 1  50 102 V12D07-075_C1      34     126       6      72      FALSE
    # # 2  14  94 V12D07-075_C1       0     120       6      70      FALSE
    # # 3  61  97 V12D07-075_C1      43     127       5      71      FALSE
    # 
    # spots$row_distance <- with(spots, pmin(abs(row - min_row), abs(row - max_row)))
    # spots$col_distance <- with(spots, pmin(abs(col - min_col), abs(col - max_col)))
    # spots$edge_distance <- with(spots, pmin(row_distance, col_distance))
    # spe$edge_spots <- factor(spots$edge_spots, levels = c("TRUE", "FALSE"))
    # spe$edge_distance <- spots$edge_distance
    # 
    # spe$scran_low_lib_size_edge <-
    #     factor(
    #         qcfilter$low_lib_size &
    #             spots$edge_distance < 1,
    #         levels = c("TRUE", "FALSE")
    #     )

    return(spe)
}

spe <- metrics_qc(spe)
colnames(colData(spe))


lobstr::obj_size(spe)
# 1.92 GB

## Save object with metrics_qc()
# saveRDS(spe, file.path(dir_rdata, "spe_qc.rds"))

## Locate low library size spots on the edge
addmargins(table("Low_libsize_edge" = spe$scran_low_lib_size_edge))
# Low_libsize_edge
# TRUE FALSE   Sum
# 31 16897 16928

# ==============================================================================
## plot low library size spots in-tissues & in the edge
## egde_distance sample plots
# edge_distance https://github.com/LieberInstitute/Visium_SPG_AD/blob/master/plots/07_spot_qc/egde_distance_wholegenome.pdf
# https://github.com/LieberInstitute/Visium_SPG_AD/blob/master/plots/07_spot_qc/scran_targeted_low_lib_size_vs_edge_distance.pdf
# https://github.com/LieberInstitute/Visium_SPG_AD/blob/master/plots/07_spot_qc/scran_targeted_low_lib_size.pdf


## -----------------------------
## scran in-tissue metrics


## low library size in-tissue edge distance

vis_grid_gene(
    spe = spe,
    geneid = "edge_distance",
    assayname = "counts",
    sample_order = lst_order,
    height = var_height,
    width = var_width,
    point_size = var_point_size,
    # return_plots = TRUE,
    pdf = here(dir_plots, "in_tissue_egde_distance.pdf"),
    spatial = FALSE,
    minCount = -1,
    cont_colors = viridisLite::viridis(21, direction = -1)
)



print("Ploting in tissue scran metrics")

lst_in_scran_counts <- c(
    scran_low_lib_size = "in_tissue_scran_low_lib_size.pdf",
    scran_low_lib_size_edge = "in_tissue_scran_low_lib_size_edge.pdf",
    scran_high_subsets_Mito_percent = "in_tissue_scran_high_Mito_percent.pdf"
)
lst_size_spot <- c((var_point_size + 1), (var_point_size + 1), (var_point_size + 1))
lst_scran_vars <- list((names(lst_in_scran_counts)), (lst_in_scran_counts), lst_size_spot)

plt_scran_func <- function(idvar, pdf_name, spot_s) {
    vis_grid_clus(
        spe = spe,
        clustervar = idvar,
        sample_order = lst_order,
        height = var_height, # 8
        width = var_width, # 9
        point_size = spot_s,
        pdf = here(dir_plots, pdf_name),
        sort_clust = FALSE,
        colors = c("TRUE" = "blue", "FALSE" = "grey90")
    )
}

pmap(lst_scran_vars, plt_scran_func)


## low library size and chrM ratios

low_library <- map(unique(spe$sample_id), ~ summary(spe$scran_low_lib_size[spe$sample_id == .x]))
print(unlist(low_library))
low_library_T <- sum(as.numeric(sapply(low_library, "[[", 1)))

message(low_library_T, " spots detected with scran_low_lib_size ")

low_library_edge <- map(unique(spe$sample_id), ~ summary(spe$scran_low_lib_size_edge[spe$sample_id == .x]))
low_library_edge_F <- sum(as.numeric(sapply(low_library_edge, "[[", 2)))
tmp <- map(seq_along(low_library_edge), ~ as.integer(low_library_edge[[.x]][3])) |> unlist() 
low_library_edge_T <-  sum(replace(tmp, is.na(tmp), 0))

message(low_library_edge_T, " spots detected with scran_low_lib_size at edge")

## Get summary for chrM ratio versus high_subsets_Mito_percent detected by scran for reference

map(unique(spe$sample_id), ~ summary(spe$expr_chrM[spe$sample_id == .x]))
map(unique(spe$sample_id), ~ sum(spe$expr_chrM[spe$sample_id == .x]))

map(unique(spe$sample_id), ~ summary(spe$expr_chrM_ratio[spe$sample_id == .x] * 100))
map(unique(spe$sample_id), ~ summary(spe$scran_high_subsets_Mito_percent[spe$sample_id == .x]))
map(unique(spe$sample_id), ~ head(spe$scran_high_subsets_Mito_percent[spe$sample_id == .x])) # boolean
# map(unique(spe$sample_id), ~ (spe$expr_chrM_ratio * 100))
# colnames(colData(spe))


# ==============================================================================
## Drop spots with a low library size that are on the edge

message("Initial number of spots:", dim(spe)[2], "\n")

spe <- spe[, spe$scran_low_lib_size_edge == "FALSE"]
message(
    "Number of spots after removed low library size spots on the tissue edge:",
    dim(spe)[2],
    "\n"
)


lobstr::obj_size(spe)
# 1.16 GB

## Second round to remove any remaining empty spots and/or genes with zero counts

spe <- spe[
    rowSums(assays(spe)$counts) > 0,
    (colSums(assays(spe)$counts) > 0) & spe$in_tissue
]
# spe1 <- spe$in_tissue[ (rowSums(assays(spe)$counts) > 0), (colSums(assays(spe)$counts) > 0) ]

message("Number of spots after preliminary QC:", dim(spe)[2], "\n")
# Number of spots after preliminary QC: 16897

# Note: saved here the object to use in the app

saveRDS(spe, file.path(dir_rdata, "spe_qc_low_lib_edge.rds"))


# ==============================================================================

## Additional QC. Drop spots with high chrM percentage marked as outlier by scran

# spe <- spe[, spe$scran_high_subsets_Mito_percent == "FALSE"]
# message(
#     "Number of spots after removed high chrM percentage spots:",
#     dim(spe)[2],
#     "\n"
# )


## Save object with metrics_qc()

# saveRDS(spe, file.path(dir_rdata, "spe_qc_low_lib_edge_HighM.rds"))



# ==============================================================================

# library("slurmjobs")
# job_single(
#   name = "02_scran_exploratory_QCs", memory = "50G", cores = 2, create_shell = TRUE
# )


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# > Sys.time()
# [1] "2024-04-16 16:06:31 EDT"
# > proc.time()
# user   system  elapsed
# 69.713    2.288 3585.754
# > options(width = 120)
# > session_info()
# ─ Session info ──────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.2 Patched (2024-02-08 r85876)
# os       Rocky Linux 9.2 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2024-04-16
# pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3.x/bin/pandoc
#
# ─ Packages ──────────────────────────────────────────────────────────────────────────────────────────────
# package                * version     date (UTC) lib source
# abind                    1.4-5       2016-07-21 [2] CRAN (R 4.3.2)
# AnnotationDbi            1.64.1      2023-11-03 [2] Bioconductor
# AnnotationHub            3.10.0      2023-10-24 [2] Bioconductor
# attempt                  0.3.1       2020-05-03 [2] CRAN (R 4.3.2)
# beachmat                 2.18.0      2023-10-24 [2] Bioconductor
# beeswarm                 0.4.0       2021-06-01 [2] CRAN (R 4.3.2)
# benchmarkme              1.0.8       2022-06-12 [2] CRAN (R 4.3.2)
# benchmarkmeData          1.0.4       2020-04-23 [2] CRAN (R 4.3.2)
# Biobase                * 2.62.0      2023-10-24 [2] Bioconductor
# BiocFileCache            2.10.1      2023-10-26 [2] Bioconductor
# BiocGenerics           * 0.48.1      2023-11-01 [2] Bioconductor
# BiocIO                   1.12.0      2023-10-24 [2] Bioconductor
# BiocManager              1.30.22     2023-08-08 [2] CRAN (R 4.3.2)
# BiocNeighbors            1.20.2      2024-01-07 [2] Bioconductor 3.18 (R 4.3.2)
# BiocParallel             1.36.0      2023-10-24 [2] Bioconductor
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
# dplyr                  * 1.1.4       2023-11-17 [2] CRAN (R 4.3.2)
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
# forcats                * 1.0.0       2023-01-29 [2] CRAN (R 4.3.2)
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
# here                   * 1.0.1       2020-12-13 [2] CRAN (R 4.3.2)
# hms                      1.1.3       2023-03-21 [2] CRAN (R 4.3.2)
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
# lobstr                 * 1.1.2       2022-06-22 [2] CRAN (R 4.3.2)
# locfit                   1.5-9.8     2023-06-11 [2] CRAN (R 4.3.2)
# lubridate              * 1.9.3       2023-09-27 [2] CRAN (R 4.3.2)
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
# prettyunits              1.2.0       2023-09-24 [2] CRAN (R 4.3.2)
# promises                 1.2.1       2023-08-10 [2] CRAN (R 4.3.2)
# purrr                  * 1.0.2       2023-08-10 [2] CRAN (R 4.3.2)
# R6                       2.5.1       2021-08-19 [2] CRAN (R 4.3.2)
# rappdirs                 0.3.3       2021-01-31 [2] CRAN (R 4.3.2)
# RColorBrewer             1.1-3       2022-04-03 [2] CRAN (R 4.3.2)
# Rcpp                     1.0.12      2024-01-09 [2] CRAN (R 4.3.2)
# RCurl                    1.98-1.14   2024-01-09 [2] CRAN (R 4.3.2)
# readr                  * 2.1.5       2024-01-10 [2] CRAN (R 4.3.2)
# rematch2                 2.1.2       2020-05-01 [2] CRAN (R 4.3.2)
# restfulr                 0.0.15      2022-06-16 [2] CRAN (R 4.3.2)
# rjson                    0.2.21      2022-01-09 [2] CRAN (R 4.3.2)
# rlang                    1.1.3       2024-01-10 [2] CRAN (R 4.3.2)
# rprojroot                2.0.4       2023-11-05 [2] CRAN (R 4.3.2)
# Rsamtools                2.18.0      2023-10-24 [2] Bioconductor
# RSQLite                  2.3.5       2024-01-21 [2] CRAN (R 4.3.2)
# rsvd                     1.0.5       2021-04-16 [2] CRAN (R 4.3.2)
# rtracklayer              1.62.0      2023-10-24 [2] Bioconductor
# S4Arrays                 1.2.0       2023-10-24 [2] Bioconductor
# S4Vectors              * 0.40.2      2023-11-23 [2] Bioconductor 3.18 (R 4.3.2)
# sass                     0.4.8       2023-12-06 [2] CRAN (R 4.3.2)
# ScaledMatrix             1.10.0      2023-10-24 [2] Bioconductor
# scales                   1.3.0       2023-11-28 [2] CRAN (R 4.3.2)
# scater                   1.30.1      2023-11-16 [2] Bioconductor
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
# spatialLIBD            * 1.14.1      2023-11-30 [2] Bioconductor 3.18 (R 4.3.2)
# statmod                  1.5.0       2023-01-06 [2] CRAN (R 4.3.2)
# stringi                  1.8.3       2023-12-11 [2] CRAN (R 4.3.2)
# stringr                * 1.5.1       2023-11-14 [2] CRAN (R 4.3.2)
# SummarizedExperiment   * 1.32.0      2023-10-24 [2] Bioconductor
# tibble                 * 3.2.1       2023-03-20 [2] CRAN (R 4.3.2)
# tidyr                  * 1.3.1       2024-01-24 [2] CRAN (R 4.3.2)
# tidyselect               1.2.0       2022-10-10 [2] CRAN (R 4.3.2)
# tidyverse              * 2.0.0       2023-02-22 [2] CRAN (R 4.3.2)
# timechange               0.3.0       2024-01-18 [2] CRAN (R 4.3.2)
# tzdb                     0.4.0       2023-05-12 [2] CRAN (R 4.3.2)
# utf8                     1.2.4       2023-10-22 [2] CRAN (R 4.3.2)
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
#
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────

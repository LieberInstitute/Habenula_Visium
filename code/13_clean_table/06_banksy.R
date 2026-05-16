library(here)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(sessioninfo)
library(Banksy)
library(cowplot)
library(scater)

lambda = 0.2
res_neat = "res1_8"

spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'banksy',
    'spe_banksy.rds'
)
spe_orig_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)

# read cluster file run before
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'banksy',
    sprintf('leiden_%s.csv', res_neat)
)

plot_dir = here(
    'plots', '09_HD_cell_level', 'no_secondary', 'banksy',
    sprintf('leiden_%s', res_neat)
)

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# helper functions for subtype-based shades
# -----------------------------
mix_color <- function(col1, col2, p = 0.5) {
    rgb1 <- col2rgb(col1) / 255
    rgb2 <- col2rgb(col2) / 255
    rgb_new <- (1 - p) * rgb1 + p * rgb2
    rgb(rgb_new[1], rgb_new[2], rgb_new[3])
}

shift_color <- function(base_col, amount) {
    if (amount < 0) {
        mix_color(base_col, "white", -amount)
    } else if (amount > 0) {
        mix_color(base_col, "black", amount)
    } else {
        base_col
    }
}

make_subtype_shades <- function(base_col, n) {
    if (n == 1) return(base_col)
    amounts <- seq(-0.45, 0.25, length.out = n)
    vapply(amounts, function(a) shift_color(base_col, a), character(1))
}

# -----------------------------
# fine subtype main colors
# -----------------------------
fine_colors = c(
    OPC = '#d3c871',
    Oligo = '#4d5802',
    Microglia = '#222222',
    Astrocyte = '#8d363c',
    Endo = '#ee6c14',
    MHb.1 = '#FF00FF',
    MHb.2 = '#FAA0A0',
    LHb.2.7 = '#00a900',
    LHb.1.3.4 = '#004F2D',
    LHb.4 = '#84DCC6',
    Excit.Thal = '#9e4ad1'
)

cluster_anno <- read.csv(
    "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/09_HD_cell_level/no_secondary/registration_banksy/cluster_annotation.csv"
)

fine_lookup = setNames(
    cluster_anno$fine_cell_type,
    as.character(cluster_anno$cluster)
)

# -----------------------------
# read data
# -----------------------------
spe = readRDS(spe_path)
cluster_df = read.csv(cluster_path)

# -----------------------------
# detect cluster column
# -----------------------------
if ("cluster" %in% colnames(cluster_df)) {
    cluster_col = "cluster"
} else if ("banksy" %in% colnames(cluster_df)) {
    cluster_col = "banksy"
} else {
    stop("Cannot find cluster column in cluster_path.")
}

if (!"key" %in% colnames(cluster_df)) {
    stop("cluster file must contain a 'key' column to match spe$key.")
}

cluster_df = cluster_df %>%
    mutate(
        key = as.character(key),
        cluster = as.character(.data[[cluster_col]])
    ) %>%
    select(key, cluster)

# check match
if (!all(spe$key %in% cluster_df$key)) {
    missing_keys = setdiff(spe$key, cluster_df$key)
    stop("Some spe$key values are missing in cluster file. Example: ",
         paste(head(missing_keys), collapse = ", "))
}

# match cluster back to spe
cluster_vec = cluster_df$cluster[match(spe$key, cluster_df$key)]
spe$cluster_plot = factor(cluster_vec, levels = sort(unique(cluster_vec)))

# -----------------------------
# detect UMAP name
# -----------------------------
rd_name = reducedDimNames(spe)[
    grep(sprintf('^UMAP.*lam%s', lambda), reducedDimNames(spe))
]

if (length(rd_name) != 1) {
    stop("Expected exactly one rd_name, found: ", paste(rd_name, collapse = ", "))
}

# -----------------------------
# fix spatial coordinates
# -----------------------------
spe_orig = readRDS(spe_orig_path)
spatialCoords(spe) = spatialCoords(spe_orig[, spe$key])
rm(spe_orig)
gc()

# -----------------------------
# build cluster palette from subtype families
# -----------------------------
cluster_ids = levels(spe$cluster_plot)
cluster_subtype = fine_lookup[cluster_ids]

if (any(is.na(cluster_subtype))) {
    stop(
        "Some cluster IDs do not have matched fine_cell_type in cluster_anno: ",
        paste(cluster_ids[is.na(cluster_subtype)], collapse = ", ")
    )
}

if (any(is.na(fine_colors[cluster_subtype]))) {
    stop("Some fine_cell_type values do not have matched colors in fine_colors.")
}

cluster_palette = setNames(rep(NA_character_, length(cluster_ids)), cluster_ids)

for (ct in unique(cluster_subtype)) {
    cl_ids = names(cluster_subtype)[cluster_subtype == ct]
    cl_ids_num = suppressWarnings(as.numeric(cl_ids))

    if (!any(is.na(cl_ids_num))) {
        cl_ids = as.character(sort(cl_ids_num))
    } else {
        cl_ids = sort(cl_ids)
    }

    sub_cols = make_subtype_shades(
        base_col = fine_colors[ct],
        n = length(cl_ids)
    )

    cluster_palette[cl_ids] = sub_cols
}

# numeric legend order
cluster_ids = unique(as.character(cluster_vec))
cluster_ids = as.character(sort(as.numeric(cluster_ids)))

spe$cluster_plot = factor(
    as.character(cluster_vec),
    levels = cluster_ids
)

cluster_palette = cluster_palette[cluster_ids]

# -----------------------------
# plot clusters and UMAP
# -----------------------------
for (sample_id in unique(spe$sample_id)) {
    if (length(unique(spe$cluster_plot)) <= 36) {
        p = vis_clus(
                spe,
                sampleid = sample_id,
                clustervar = "cluster_plot",
                is_stitched = TRUE,
                point_size = 20,
                spatial = FALSE
            ) +
            scale_fill_manual(
                values = cluster_palette,
                breaks = cluster_ids,
                drop = FALSE
            ) +
            guides(fill = guide_legend(override.aes = list(size = 8)))

        png(
            file.path(plot_dir, sprintf('clusters_%s_new.png', sample_id)),
            width = 1500, height = 1500
        )
        print(p)
        dev.off()
    }

    p = plotReducedDim(
            spe[, spe$sample_id == sample_id],
            dimred = rd_name,
            point_size = 0.6,
            colour_by = "cluster_plot"
        ) +
        scale_colour_manual(
            values = cluster_palette,
            breaks = cluster_ids,
            drop = FALSE
        ) +
        theme_bw(base_size = 20) +
        guides(color = guide_legend(override.aes = list(size = 4)))

    png(
        file.path(plot_dir, sprintf('UMAP_%s_new.png', sample_id)),
        width = 1500, height = 1500
    )
    print(p)
    dev.off()
}

session_info()
library(here)
library(tidyverse)
library(sessioninfo)
library(spatialLIBD)
library(duckplyr)

cor_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
    'registration', 'cor_rds', 'cleaning_y', 'cor_vs_snRNAseq_fine.rds'
)
cluster_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
    'k14_cluster_coords.csv.gz'
)
plot_dir = here(
    'plots', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
    'cleany', 'k_14_manual'
)
cor_index = 12
manual_anno = c(
    '2' = 'MHb.1.2',
    '9' = 'Ambig_Microglia_1',
    '5' = 'Ambig_Microglia_2'
)
cell_type_colors = c(
    OPC = '#d3c871',
    Oligo = '#4d5802',
    Microglia = '#222222',
    Astrocyte = '#8d363c',
    Endo = '#ee6c14',
    MHb.1 = '#FF00FF',
    MHb.1.2 = '#D9678B',
    MHb.2 = '#FAA0A0',
    LHb.2.7 = '#00a900',
    LHb.1.3.4 = '#004F2D',
    LHb.4 = '#84DCC6',
    Excit.Thal = '#9e4ad1',
    Ambig_LHb.7 = '#494949',
    Ambig_Microglia_1 = '#616161',
    Ambig_Microglia_2 = '#989898',
    Ambig_Oligo = '#CCCCCC'
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 1))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(plot_dir, showWarnings = FALSE)

#   Using the spatial registration for k = 14 (optimal result), annotate
#   FICTURE clusters with cell types. Explicitly label ambiguous clusters since
#   most have very dirty matches to cell types
anno_df = readRDS(cor_path)[[cor_index]] |>
    annotate_registered_clusters(cutoff_merge_ratio = 0.1) |>
    as_tibble() |>
    mutate(
        cell_type = case_when(
            cluster %in% names(manual_anno) ~ manual_anno[cluster],
            grepl('\\*$', layer_label) ~ str_replace(
                layer_label, '(.+)\\*$', 'Ambig_\\1'
            ),
            TRUE ~ layer_label
        ),
        cluster = as.integer(cluster)
    ) |>
    select(cluster, cell_type) |>
    arrange(cluster)

#   Read in pre-computed FICTURE cluster calls with the associated spatial
#   coordinates
cluster_df = read_csv_duckdb(cluster_path, prudence = 'stingy') |>
    collect() |>
    mutate(
        cell_type = factor(
            anno_df$cell_type[match(as.integer(FICTURE_k14), anno_df$cluster)],
            levels = names(cell_type_colors)
        )
    )
stopifnot(!any(is.na(cluster_df$cell_type)))

#   Plot each donor separately
for (sample_id in unique(cluster_df$sample_id)) {
    p = cluster_df |>
        filter(sample_id == !!sample_id) |>
        ggplot(aes(x = X, y = 0 - Y, color = cell_type)) +
            geom_point(size = 0.2, shape = 15) +
            scale_color_manual(values = cell_type_colors) +
            coord_fixed() +
            theme_void(base_size = 100) +
            labs(color = 'Cell Type') +
            guides(color = guide_legend(override.aes = list(size = 10)))
    pdf(file.path(plot_dir, paste0(sample_id, '.pdf')), width = 80, height = 80)
    print(p)
    dev.off()
}

session_info()

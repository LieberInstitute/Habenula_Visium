library(tidyverse)
library(here)
library(sessioninfo)
library(duckplyr)

banksy_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'banksy',
    'leiden_res1_8.csv'
)
ficture_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'ficture_plotting', 'extracellular.parquet'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'no_secondary', 'misc_paper_figs'
)
donor_colors = c(
    Br9090 = '#5465FF',
    Br8433 = '#474B24',
    Br8667 = '#CEF9F2',
    Br3942 = '#3D0B37',
    Br9902 = '#D6CA98'
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

composition_barplot = function(this_df, x_lab, y_lab, plot_path) {
    p = this_df |>
        count(cluster, donor) |>
        group_by(cluster) |>
        mutate(prop = n / sum(n)) |>
        ungroup() |>
        mutate(cluster = factor(cluster, levels = sort(unique(cluster)))) |>
        ggplot(aes(x = cluster, y = prop, fill = donor)) +
            geom_col() +
            scale_fill_manual(values = donor_colors) +
            scale_y_continuous(labels = scales::percent_format()) +
            labs(x = x_lab, y = y_lab, fill = "Donor") +
            theme_bw(base_size = 16)

    pdf(plot_path, width = 10, height = 4)
    print(p)
    dev.off()
}

ficture_df = read_parquet_duckdb(ficture_path, prudence = 'stingy') |>
    dplyr::rename(cluster = k17) |>
    filter(!is.na(cluster)) |>
    select(bin_key, cluster) |>
    collect() |>
    mutate(
        donor = factor(
            paste0('Br', str_extract(bin_key, '[0-9]{4}$')),
            levels = names(donor_colors)
        )
    ) |>
    select(donor, cluster)
stopifnot(!any(is.na(ficture_df$donor)))

banksy_df = read_csv_duckdb(banksy_path, prudence = 'stingy') |>
    collect() |>
    mutate(
        donor = factor(
            paste0('Br', str_extract(key, '[0-9]{4}$')),
            levels = names(donor_colors)
        )
    ) |>
    dplyr::rename(cluster = banksy) |>
    select(donor, cluster)

composition_barplot(
    ficture_df, "Extracellular FICTURE cluster", "Proportion of bins",
    file.path(plot_dir, 'ficture_donor_composition.pdf')
)
composition_barplot(
    banksy_df, "Banksy cluster", "Proportion of cells",
    file.path(plot_dir, 'banksy_donor_composition.pdf')
)

session_info()

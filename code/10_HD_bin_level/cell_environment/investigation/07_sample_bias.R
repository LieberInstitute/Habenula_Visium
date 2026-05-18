#   How sample-biased are the clustering results for extracellular and all bins?
#   I interactively tried a bunch of k values and surprsingly, extracellular
#   results are not systematically more biased than all-bin results despite
#   sample bias being much stronger at the gene-expression level. Also,
#   subsetting all-bin results to extracellular bins doesn't make bias worse
#   (tested that interactively)

library(here)
library(sessioninfo)
library(tidyverse)
library(duckplyr)

k = 10
extra_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'ficture_plotting', 'extracellular.parquet'
)

all_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
    'bin_level_clusters_batch.parquet'
)
plot_dir = here(
    'plots', '10_HD_bin_level', 'no_secondary', 'cell_environment', 
    'investigation'
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

extra_df = read_parquet_duckdb(extra_path, prudence = 'lavish') |>
    mutate(
        sample_id = stringr::str_extract(bin_key, "H1[^_]*_[A-Z][0-9]_[0-9]{4}")
    ) |>
    dplyr::rename(ficture_cluster = paste0('k', k)) |>
    filter(!is.na(ficture_cluster)) |>
    select(sample_id, ficture_cluster)

all_df = read_parquet_duckdb(all_path, prudence = 'lavish') |>
    dplyr::rename(ficture_cluster = paste0('FICTURE_k', k)) |>
    filter(!is.na(ficture_cluster), ficture_cluster != "NA") |>
    select(sample_id, ficture_cluster)

# ## Compositional barplots: cluster x sample proportion, faceted by dataset
plot_df = bind_rows(
        extra_df |>
            count(sample_id, cluster = factor(ficture_cluster)) |>
            mutate(prop = n / sum(n), .by = cluster) |>
            mutate(dataset = "extracellular"),
        all_df |>
            filter(!is.na(ficture_cluster)) |>
            count(sample_id, cluster = factor(ficture_cluster)) |>
            mutate(prop = n / sum(n), .by = cluster) |>
            mutate(dataset = "all bins")
    ) |>
    collect() |>
    mutate(dataset = factor(dataset, levels = c("extracellular", "all bins")))

p = ggplot(plot_df, aes(x = cluster, y = prop, fill = sample_id)) +
    geom_col() +
    facet_wrap(~dataset, nrow = 2) +
    scale_y_continuous(labels = scales::percent) +
    labs(
        x = paste0("FICTURE cluster (k=", k, ")"),
        y = "Percentage",
        fill = "Sample",
        title = "Sample composition by cluster"
    ) +
    theme_bw(base_size = 15)

pdf(file.path(plot_dir, sprintf('sample_composition_k%d.pdf', k)))
print(p)
dev.off()

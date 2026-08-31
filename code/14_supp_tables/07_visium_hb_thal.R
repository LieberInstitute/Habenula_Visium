library(here)
library(tidyverse)
library(sessioninfo)
library(SpatialExperiment)
library(spatialLIBD)

spe_in_path = here("processed-data", "04_harmony_BayesSpace", "spe_harmony.rds")
cluster_path = here(
    "processed-data", "04_harmony_BayesSpace", "clusters_BayesSpace",
    "BayesSpace_harmony_k28", "clusters.csv"
)
plot_dir = here("plots", "14_supp_tables")
hb_clusters = c(
    '1' = '#AC0561',
    '5' = '#4357AD',
    '10' = '#594236',
    '11' = '#058837',
    '20' = '#DF9526',
    '27' = '#48A9A6',
    'Other' = '#9A9A9A'
)
thal_clusters = c(
    '17' = '#1F487E',
    '21' = '#F78E69',
    'Other' = '#9A9A9A'
)
#   Selected for good counts of both hb and thal clusters
sample_id = 'V13B23-280_C1'

spe = readRDS(spe_in_path)
spe$bayesspace = tibble(key = spe$key) |>
    left_join(read_csv(cluster_path), by = "key") |>
    pull(cluster)

spe$hb = factor(
    ifelse(spe$bayesspace %in% names(hb_clusters), spe$bayesspace, "Other"),
    levels = names(hb_clusters)
)
spe$thal = factor(
    ifelse(spe$bayesspace %in% names(thal_clusters), spe$bayesspace, "Other"),
    levels = names(thal_clusters)
)

for (region in c('hb', 'thal')) {
    #   Bug where first plot of the session looks different than all others
    for (i in seq_len(2)) {
        p = vis_clus(
            spe, sampleid = sample_id, clustervar = region,
            colors = get(paste0(region, '_clusters')), spatial = FALSE
        )
    }
    pdf(file.path(plot_dir, sprintf('visium_%s.pdf', region)))
    print(p)
    dev.off()
}

session_info()

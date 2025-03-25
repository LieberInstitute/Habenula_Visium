#   On the SPE object pseudobulked by Banksy cluster at Leiden resolution 1,
#   perform PCA to see if points separate by probe set lot group, which would
#   suggest a type of batch effect described by 10X

library(sessioninfo)
library(here)
library(SpatialExperiment)
library(scran)
library(scater)
library(BiocSingular)
library(tidyverse)

spe_path = here(
    'processed-data', '09_HD_cell_level', 'registration_banksy',
    'pseudobulk_spe', '1.rds'
)
plot_path = here(
    'plots', '09_HD_cell_level', 'registration_banksy', 'lot_PCA.pdf'
)
num_pcs = 2

set.seed(0)

spe = readRDS(spe_path)

#   Perform PCA. Use IrlbaParam() for speed and memory,
#   inspired by https://pachterlab.github.io/voyager/articles/vig6_merfish.html#pca-for-larger-datasets
message(Sys.time(), " | Running PCA...")
spe = runPCA(spe, ncomponents = num_pcs, BSPARAM = IrlbaParam())

pc_df = tibble(
    PC1 = unname(reducedDims(spe)$PCA[, 'PC1']),
    PC2 = unname(reducedDims(spe)$PCA[, 'PC2']),
    lot = ifelse(
        spe$sample_id == 'H1-W369TJK_D1_9090',
        'Group 1', 'Group 2'
    )
)

p = ggplot(pc_df, aes(x = PC1, y = PC2, color = lot)) +
    geom_point() +
    theme_bw(base_size = 20)
pdf(plot_path)
print(p)
dev.off()

session_info()

library(here)
library(spatialLIBD)
library(sessioninfo)
library(tidyverse)
library(jaffelab)
library(Matrix)

k = 10
spe_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'registration', 'pseudobulk_spe', sprintf('%s.rds', k)
)
svg_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'nnSVG_out',
    'merged_SVGs.txt'
)

spe = readRDS(spe_path)

#   Regress out sample ID
mod = with(colData(spe), model.matrix(~ sample_id))
assays(spe)$logcounts = cleaningY(assays(spe)$logcounts, mod, P = 1)

#   Subset to SVGs for speed. We aren't guaranteed SVGs provide good signal in
#   the extracellular bins, but we need to reduce runtime somehow
svg = readLines(svg_path)
stopifnot(all(svg %in% rownames(spe)))

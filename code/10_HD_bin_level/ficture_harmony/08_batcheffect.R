library(jaffelab)
library(here)
library(tidyverse)
library(scran)
library(spatialLIBD)
library(sessioninfo)
library(HDF5Array)

spe <- readRDS("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/ficture_harmony/spe_raw.rds")
spe_raw_dir = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/ficture_harmony/spe"
message(Sys.time(), " | Saving raw SPE")
spe <- saveHDF5SummarizedExperiment(
    spe, dir = spe_raw_dir, replace = TRUE, as.sparse = TRUE
)
message(Sys.time(), " | Done saving raw SPE")

## Define a full model

# spe_in_dir = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/ficture_harmony/spe/"
# spe = loadHDF5SummarizedExperiment(spe_in_dir)

# mod <- with(colData(spe), model.matrix(~ sample_id))
# y_clean_p2 <- cleaningY(assays(spe)$counts, mod, P=2)


library("spatialLIBD")
library("here")


## Removes not required assays to deploy the shiny app

here::here("code", "03_spatialLIBD_app_pseudobulk")

path_in <- here("processed-data", "04_harmony_BayesSpace", "spe_harmony.rds")
# path_out <- here("processed-data", "04_harmony_BayesSpace", "spe_harmony_shiny.rds")
path_out <- here("processed-data", "05_brain_area_differential_expression", "spe_pseudobulk_shiny.rds")


## symbolic link to point the spe.rds object to clean
spe <- readRDS(path_in) # spe with harmony and BayesSpace

# lobstr::obj_size(spe)
# 10.89 GB with 5 capture areas
# lobstr::obj_size(reducedDim(spe, "TSNE_perplexity05.harmony_subject_no_lambda"))
# lobstr::obj_size(assay(spe, "counts"))

# Removing not required assays to reduce size of instance
assay(spe, "counts") <- NULL
assay(spe, "binomial_deviance_residuals") <- NULL
lobstr::obj_size(spe)
# 1.37 GB

saveRDS(spe, path_out)

## Set up soft links if needed
withr::with_dir(
    here("code", "03_spatialLIBD_app_pseudobulk"),
    system("ln -s ../../processed-data/04_harmony_BayesSpace/spe_pseudobulk_shiny.rds spe_pseudobulk_shiny.rds")
)


library("spatialLIBD")
library("here")


## This removes not required assays to deploy the shiny app

here("code", "03_spatialLIBD_app")
path_in <- path_out <- here("processed-data", "04_harmony_BayesSpace", "spe_harmony.rds")
path_out <- here("processed-data", "04_harmony_BayesSpace", "spe_harmony_shiny.rds")

## symbolic link to point the spe.rds object to clean
spe <- readRDS(path_in) # spe with harmony and BayesSpace

# lobstr::obj_size(spe)
# 10.89 GB with 5 capture areas
# lobstr::obj_size(reducedDim(spe, "TSNE_perplexity05.harmony_subject_no_lambda"))
# lobstr::obj_size(assay(spe, "counts"))

# Removing not required assays to reduce size of instance
assay(spe, "counts") <- NULL
assay(spe, "binomial_deviance_residuals") <- NULL
# lobstr::obj_size(spe)
# 4.03 GB

saveRDS(spe, path_out)


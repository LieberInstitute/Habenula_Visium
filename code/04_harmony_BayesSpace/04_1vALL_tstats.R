library("here")
library("spatialLIBD")
stopifnot(packageVersion("spatialLIBD") >= "1.9.18")
library("sessioninfo")
library("scater")

## output directory
dir_rdata <- here::here("processed-data", "11_grey_matter_only", opt$spetype)
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dir_rdata)) ## Check that it was created successfully

## load spe data
spe <-
    readRDS(
        here::here(
            "processed-data",
            "08_harmony_BayesSpace",
            opt$spetype,
            paste0("spe_harmony_", opt$spetype, ".rds")
        )
    )

spe <- cluster_import(
    spe,
    cluster_dir = here::here(
        "processed-data",
        "08_harmony_BayesSpace",
        opt$spetype,
        "clusters_BayesSpace"
    ),
    prefix = ""
)

spe <- cluster_import(
    spe,
    cluster_dir = here::here(
        "processed-data",
        "09_pathology_vs_BayesSpace",
        "pathology_levels"
    ),
    prefix = ""
)

library("rsconnect")
library("here")

## Or you can go to your shinyapps.io account and copy this
## Here we do this to keep our information hidden.
# load(here("code", "03_spatialLIBD_app", ".deploy_info.Rdata"), verbose = TRUE)
# rsconnect::setAccountInfo(
#     name = deploy_info$name,
#     token = deploy_info$token,
#     secret = deploy_info$secret
# )

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Set BayesSpace k selection
BayesSpace_k <- 24


## Build names
appDir_name <-  here("code", "03_spatialLIBD_app_pseudobulk")
sce_pseudo_name <- paste0("sce_pseudo_BayesSpace_k", BayesSpace_k, ".rds")
modeling_results_name <- paste0("modeling_results_BayesSpace_k", BayesSpace_k, ".Rdata")
sig_genes_name <- paste0("sig_genes_k", BayesSpace_k, ".Rdata")
app_name <- paste0("Habenula_Visium_Sp", BayesSpace_k)

## Deploy the app, that is, upload it to shinyapps.io
rsconnect::deployApp(
    appDir = appDir_name,
    appFiles = c(
        "app.R",
        "spe_subset_for_spatialLIBD.rds",
        sce_pseudo_name,
        modeling_results_name,
        sig_gene,
        withr::with_dir(appDir_name, dir("clusters_BayesSpace", full.names = TRUE)),
        withr::with_dir(appDir_name, dir("www", full.names = TRUE))
    ),
    appName = app_name,
    account = "libd",
    server = "shinyapps.io"
)

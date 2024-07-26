## copied from https://github.com/LieberInstitute/spatialDLPFC/blob/bd93c980d7653579f81ff1c91c309cea0c7474a6/code/deploy_app_k09/subset.R

library("spatialLIBD")
library("lobstr")
library("here")
library("sessioninfo")


## Set working directory
here("code", "03_spatialLIBD_app_pseudobulk")

## Set BayesSpace k selection
BayesSpace_k <- 24

## Set up soft links if needed
model_result <- paste0('modeling_results_BayesSpace_k', BayesSpace_k, '.Rdata')
softLink_command <- paste0('ln -s ../../processed-data/05_layer_differential_expression/modeling_results_BS/', model_result, ' ', model_result)
withr::with_dir(
    here("code", "03_spatialLIBD_app_pseudobulk"),
    system(softLink_command)
)

pseudo_result <- paste0('sce_pseudo_BayesSpace_k', BayesSpace_k, '.rds')
softLink_command <- paste0('ln -s ../../processed-data/05_layer_differential_expression/', pseudo_result, ' ', pseudo_result)
withr::with_dir(
    here("code", "03_spatialLIBD_app_pseudobulk"),
    system(softLink_command)
)

sig_genes_result <- paste0('sig_genes_k', BayesSpace_k, '.Rdata')
softLink_command <- paste0('ln -s ../../processed-data/05_layer_differential_expression/', sig_genes_result, ' ', sig_genes_result)
withr::with_dir(
    here("code", "03_spatialLIBD_app_pseudobulk"),
    system(softLink_command)
)

withr::with_dir(
    here("code", "03_spatialLIBD_app_pseudobulk"),
    system("ln -s ../../processed-data/04_harmony_BayesSpace/spe_harmony_shiny.rds spe_subset_for_spatialLIBD.rds")
)

withr::with_dir(
  here("code", "03_spatialLIBD_app_pseudobulk", "www"),
  system("ln -s ../../../README.md README.md")
)


## load the pseudobulked object sce_pseudo
sce_pseudo <- readRDS(pseudo_result)
# lobstr::obj_size(sce_pseudo) # 6.21 MB

## load modeling results for k16
load(model_result, verbose = TRUE)
model_result
# lobstr::obj_size(modeling_results) # 17.94 MB

## For sig_genes_extract_all() to work https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/05_deploy_app_wholegenome/app.R#L37-L44
## Quick inspection of the BayesSpace k selection
names(colData(sce_pseudo))
length(levels(sce_pseudo$BayesSpace))

## Extract BayesSpace information to shiny spatialLIBD column
sce_pseudo$spatialLIBD <- sce_pseudo$BayesSpace

## Check that we have the right number of tests
k <- BayesSpace_k
tests <- lapply(modeling_results, function(x) {
  colnames(x)[grep("stat", colnames(x))]
})
stopifnot(length(tests$anova) == 1) ## assuming only "all"
stopifnot(length(tests$enrichment) == k)
stopifnot(length(tests$pairwise) == choose(k, 2))

## This function combines the output of sig_genes_extract() from all the layer-level (group-level) modeling results
##       and builds the data required for functions such as layer_boxplot()
sig_genes <- sig_genes_extract_all(
  n = nrow(sce_pseudo), # use `x` top ranked genes (g.e. 100) if there are memory constrains, otherwise keep all significant genes
  modeling_results = modeling_results,
  sce_layer = sce_pseudo
)
## Quick inspection
str(sig_genes)
table(sig_genes@listData$model_type)
# anova enrichment   pairwise
# 2618      62832    1445136
# table(sig_genes$test)

## Check that we have the right number of tests.
## the + 1 at the end assumes only "all"
stopifnot(length(unique(sig_genes$test)) == choose(k, 2) * 2 + k + 1)

lobstr::obj_size(sig_genes)
# 3.69 GB

## Drop parts we don't need to reduce the memory
# sig_genes
sig_genes$in_rows <- NULL
sig_genes$in_rows_top20 <- NULL
lobstr::obj_size(sig_genes)
# 181.43 MB

## Extract FDR < 5%
## From
## https://github.com/LieberInstitute/brainseq_phase2/blob/be2b7f972bb2a0ede320633bf06abe1d4ef2c067/supp_tabs/create_supp_tables.R#L173-L181
fix_csv <- function(df) {
  for (i in seq_len(ncol(df))) {
    if (any(grepl(",", df[, i]))) {
      message(paste(Sys.time(), "fixing column", colnames(df)[i]))
      df[, i] <- gsub(",", ";", df[, i])
    }
  }
  return(df)
}

z <- fix_csv(as.data.frame(subset(sig_genes, fdr < 0.05)))
## Quick inspection
colnames(z)
dim(z)  # For BayesSpace k24: [1] 212102     11
table(sig_genes$model_type)
table(z$model_type)

## Save ALL the significant DEG that passed the fdr at 5%
model_topAll_csv <- paste0("spatialHb_model_results_k", k, "_ALL.csv")
model_topAll_csv <- here("processed-data", "05_layer_differential_expression", model_topAll_csv)
write.csv(z, model_topAll_csv)

## Save a subset of the top 25 significant DEG that passed the fdr at 5%. Use in case of memory constraints
z <- subset(z, top <= 25)
## Quick inspection
dim(z)  # For BayesSpace k24: [1] 9485   11
table(z$model_type)

model_top25_csv <- paste0("spatialHb_model_results_k", k, "_FDR5perc_top25.csv")
model_top25_csv <- here("processed-data", "05_layer_differential_expression", model_top25_csv)
write.csv(z, model_top25_csv)

sig_genes_csv <- paste0("sig_genes_k" , k,".Rdata")
sig_genes_csv <- here("processed-data", "05_layer_differential_expression", sig_genes_csv)
save(sig_genes, file = sig_genes_csv)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

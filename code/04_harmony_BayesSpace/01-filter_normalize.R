library("spatialLIBD")
library("here")
library("tidyverse")
library("scran")
library("BiocParallel")
library("scater")
library("scry")
library("BiocSingular")
library("sessioninfo")
library("HDF5Array")

dir_rdata <- here("processed-data", "04_harmony_BayesSpace")
# raw_in_path <- here("processed-data", "02_build_spe", "spe.rds")
filtered_in_path <- here("processed-data", "02_build_spe", "spe_qc_low_lib_edge.rds") # This spe was processed in ~/code/02_build_spe/*
filtered_ordinary_path <- file.path(dir_rdata, "spe_filtered.rds")
filtered_hdf5_dir <- file.path(dir_rdata, "spe_filtered_hdf5")
dir_plots <- here("plots", "04_harmony_BayesSpace")
num_red_dims <- 50

num_cores <- 2 # Sys.getenv('SLURM_CPUS_ON_NODE')
set.seed(20240223)

## Create output directories
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

## load a filtered spe object
spe <- readRDS(filtered_in_path)

cat("Number of spots after removed any remaining empty spots and/or genes with zero counts:", dim(spe)[2],"\n")



################################################################################
#   Manually selection of spots to drop (issue #7)
################################################################################

#  Re-upload the spots manually annotated to resume the work as noted here https://github.com/LieberInstitute/spatialLIBD/blob/77a5303f91edb7b9ffb1ce00b4193dae5d16a8a1/R/app_server.R#L1118-L1152)

spatialLIBD_ann_file <- here("processed-data", "03_spatialLIBD_app", 
                             "spatialLIBD_ManualAnnotation_2024-05-16.csv")

# rv <- reactiveValues(ManualAnnotation = rep("NA", ncol(spe)), ContCount = data.frame(key = spe$key, COUNT = NA))
rv <- list(ManualAnnotation = rep("NA", ncol(spe)), key = spe$key)
names(rv)
head(rv$ManualAnnotation) # NA
tail(rv$key)
# key COUNT
# 16892 TTGTTAGCAAATTCGA-1_V13B23-285_D1    NA
# 16893 TTGTTCAGTGTGCTAC-1_V13B23-285_D1    NA
# 16894 TTGTTGTGTGTCAAGA-1_V13B23-285_D1    NA

unique(rv$ManualAnnotation) # NA
#unique(rv$key) # NA
length(rv$ManualAnnotation) #  16897

# read the annotation file made with the spatialLIBD shiny app  
previous_work <-
  read.csv(
    spatialLIBD_ann_file,
    header = TRUE,
    stringsAsFactors = FALSE,
    na.strings = ""
  )

head(previous_work) #list
# sample_id          spot_name     ManualAnnotation
# 1  V13B23-285_B1 AACGAAAGTCGTCCCA-1 Tissue_rolls_low_lib
# 2  V13B23-285_B1 AACGTTATCAGCACCT-1 Tissue_rolls_low_lib
# 3  V13B23-285_B1 ACATAAGTCGTGGTGA-1 Tissue_rolls_low_lib

## Update the non-NA
previous_work <-
  subset(previous_work, ManualAnnotation != "NA")

unique(previous_work['sample_id'])

# add a unique keys identificator 
previous_work$key <-
  paste0(
    previous_work$spot_name,
    "_",
    previous_work$sample_id
  )
previous_work['key']
# 1  AACGAAAGTCGTCCCA-1_V13B23-285_B1
# 2  AACGTTATCAGCACCT-1_V13B23-285_B1
# 3  ACATAAGTCGTGGTGA-1_V13B23-285_B1

# match de unique IDs and get the index row from spe
m <- match(previous_work$key, spe$key)
m
# [1]  7490  7503  7644  7667  7758  7770  7788  7792  7803  7943  7983  7984
# [13]  8033  8179  8229  8284  8402  8478  8526  8539  8568  8596  8626  8925
# ...

# set and transfer the label
#rv$ManualAnnotation[m[!is.na(m)]] <- previous_work$ManualAnnotation[!is.na(m)]
spe$ManualAnnotation[m[!is.na(m)]] <- previous_work$ManualAnnotation[!is.na(m)]
spe$key[m[!is.na(m)]] 
# [1] "AACGAAAGTCGTCCCA-1_V13B23-285_B1" "AACGTTATCAGCACCT-1_V13B23-285_B1"
# [3] "ACATAAGTCGTGGTGA-1_V13B23-285_B1" "ACCATCCGCCAACTAG-1_V13B23-285_B1"
# [5] "ACTCGATGTATTTCAT-1_V13B23-285_B1" "ACTGCTCGGAAGGATG-1_V13B23-285_B1"

#unique(rv$ManualAnnotation)
lst_manual_ann <- as.list(unique(spe$ManualAnnotation))
lst_manual_ann <- lst_manual_ann[! lst_manual_ann%in% c('NA')]

# manual_ann <- cluster_export(spe, "ManualAnnotation")  #Note.Other alternative

## Additional QC. Drop spots with manual annotations 
#spe$key[m[63]] # TTGTGAGGCATGACGC-1_V13B23-285_C1
colnames(colData(spe))


for (ann in lst_manual_ann) {
  print(paste0("Removing spots for `", ann, "` manual annotation"))
  spe <- spe[, !spe$ManualAnnotation == ann]
}
# [1] "Removing spots for `low_lib_manual` manual annotation"
# [1] "Removing spots for `high_umi_manual` manual annotation"
# [1] "Removing spots for `tissue_roll` manual annotation"
# [1] "Removing spots for `high_MTr` manual annotation"

unique(spe$ManualAnnotation)
cat("Number of spots after removed low library size spots on the tissue edge:", dim(spe)[2],"\n")
  
## Double check any remaining empty spots and/or genes with zero counts
spe <- spe[
  rowSums(assays(spe)$counts) > 0,
  (colSums(assays(spe)$counts) > 0) & spe$in_tissue
]
cat("Number of spots after removed any remaining empty spots and/or genes with zero counts:", dim(spe)[2],"\n")

## Save new spe object with spots manually annotated drop
saveRDS(spe, file.path(dir_rdata, "spe_qc_low_spatialLIBD.rds"))



################################################################################
#   Compute log-normalized counts
################################################################################



#   Filter SPE: take only spots in tissue, drop spots with 0 counts for all
#   genes, and drop genes with 0 counts in every spot
message(Sys.time(), " - Running quickCluster()")

Sys.time()
spe$scran_quick_cluster <- quickCluster(
    spe,
    BPPARAM = MulticoreParam(num_cores),
    block = spe$sample_id,
    block.BPPARAM = MulticoreParam(num_cores)
)
Sys.time()

message(Sys.time(), " - Running computeSumFactors()")
Sys.time()
spe <- computeSumFactors(spe,
    # clusters = spe$scran_quick_cluster,
    BPPARAM = MulticoreParam(num_cores)
)
Sys.time()

print("Quick cluster table:")
table(spe$scran_quick_cluster)
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 337  342  320  275  497  317  423  160  150  206  380  126  191  720  372  193 
# 17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32 
# 1157  815  367  508  909  600  151  187  743  680  139  901  612  193  663  488 
# 33   34   35   36   37 
# 261  950  475  143  804  

message(Sys.time(), " - Running checking sizeFactors()")
summary(sizeFactors(spe))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.000042  0.180453  0.619294  1.000000  1.400150 14.932511 

message(Sys.time(), " - Running logNormCounts()")
spe <- logNormCounts(spe)
# assays(2): counts logcounts

# #   Save a copy of the SPE with HDF5-backed assays, which will be important to
# #   control memory consumption later
message(Sys.time(), " - Saving HDF5-backed object to control memory later")
spe = saveHDF5SummarizedExperiment(
    spe, dir = paste0(filtered_hdf5_dir, '_temp'), replace = TRUE
)
gc()

## Save new spe object with spots manually annotated drop
saveRDS(spe, file.path(dir_rdata, "spe_qc_filtered_logcounts.rds"))



######### Move since here the code to other script: CSC


################################################################################
#   Compute PCA
################################################################################

message(Sys.time(), " - Running modelGeneVar()")
## From
## http://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html#4_variance_modelling
dec <- modelGeneVar(spe,
    block = spe$sample_id,
    BPPARAM = MulticoreParam(num_cores)
)

#plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
#curve(metadata(dec)$trend(x), col="blue", add=TRUE)

pdf(file.path(dir_plots, "scran_modelGeneVar.pdf"), useDingbats = FALSE)
mapply(function(block, blockname) {
    plot(
        block$mean,
        block$total,
        xlab = "Mean log-expression",
        ylab = "Variance",
        main = blockname
    )
    #points(metadata(block)$mean, metadata(block)$var, col="red")
    curve(metadata(block)$trend(x),
        col = "blue",
        add = TRUE
    )
}, dec$per.block, names(dec$per.block))
dev.off()

message(Sys.time(), " - Running getTopHVGs()")
top.hvgs.p1 <- getTopHVGs(dec, prop = 0.1)
print(paste("Num HVGs for top 10% prop:", length(top.hvgs.p1)))
top.hvgs.p2 <- getTopHVGs(dec, prop = 0.2)
print(paste("Num HVGs for top 20% prop:", length(top.hvgs.p2)))
top.hvgs.p5 <- getTopHVGs(dec, prop = 0.5)
print(paste("Num HVGs for top 50% prop:", length(top.hvgs.p5)))

top.hvgs.fdr5 <- getTopHVGs(dec, fdr.threshold = 0.05)
print(paste("Num HVGs at FDR = 0.05:", length(top.hvgs.fdr5)))

top.hvgs.fdr1 <- getTopHVGs(dec, fdr.threshold = 0.01)
print(paste("Num HVGs at FDR = 0.01:", length(top.hvgs.fdr1)))

save(
    top.hvgs.p1,
    top.hvgs.p2,
    top.hvgs.p5,
    top.hvgs.fdr5,
    top.hvgs.fdr1,
    file = file.path(dir_rdata, "top.hvgs.Rdata")
)

message(Sys.time(), " - Running runPCA()")
Sys.time()

# lst_top_hvgs <- list(PCA = top.hvgs.fdr5, PCA_p1 = top.hvgs.p1, PCA_p2 = top.hvgs.p2, PCA_p5 = top.hvgs.p5)
# for (x in 1:length(lst_top_hvgs)) {
#   message(' Processing ', names(lst_top_hvgs[x]))
#   spe <- runPCA(spe,
#       subset_row = hvgs,
#       ncomponents = num_red_dims,
#       name = names(lst_top_hvgs[x])
#   )
# }

spe <-
    runPCA(spe,
        subset_row = top.hvgs.fdr5,
        ncomponents = num_red_dims,
        name = "PCA"
    )
spe <-
    runPCA(spe,
        subset_row = top.hvgs.p1,
        ncomponents = num_red_dims,
        name = "PCA_p1"
    )
spe <-
    runPCA(spe,
        subset_row = top.hvgs.p2,
        ncomponents = num_red_dims,
        name = "PCA_p2"
    )
spe <-
    runPCA(spe,
        subset_row = top.hvgs.p5,
        ncomponents = num_red_dims,
        name = "PCA_p5"
    )
Sys.time()

#   Plot variance explained

# lst_PCA_elbow <- list(PCA_p1 = 'pca_elbow_p1.pdf', PCA_p2 = 'pca_elbow_p2.pdf', PCA_p5 = 'pca_elbow_p5.pdf')
# plt_elbow <- function(dim_n, plt_name) {
#   percent.var <- attr(reducedDim(spe, dim_n), "percentVar")
#   pdf(file.path(dir_plots, plt_name), useDingbats = FALSE)
#   plot(percent.var, xlab = gsub("^PCA_", "PC_", dim_n), ylab = "Variance explained (%)")
#   dev.off()
# }
# 
# map2(names(lst_PCA_elbow), lst_PCA_elbow, ~ plt_elbow(.x, .y))

percent.var <- attr(reducedDim(spe, "PCA_p1"), "percentVar")
pdf(file.path(dir_plots, "pca_elbow_p1.pdf"), useDingbats = FALSE)
plot(percent.var, xlab = "PC_p1", ylab = "Variance explained (%)")
dev.off()
#
percent.var <- attr(reducedDim(spe, "PCA_p2"), "percentVar")
pdf(file.path(dir_plots, "pca_elbow_p2.pdf"), useDingbats = FALSE)
plot(percent.var, xlab = "PC_p2", ylab = "Variance explained (%)")
dev.off()
#
percent.var <- attr(reducedDim(spe, "PCA_p5"), "percentVar")
pdf(file.path(dir_plots, "pca_elbow_p5.pdf"), useDingbats = FALSE)
plot(percent.var, xlab = "PC_p5", ylab = "Variance explained (%)")
dev.off()


################################################################################
#   Compute GLM-PCA
################################################################################

message(Sys.time(), " - Running devianceFeatureSelection()")
spe <- devianceFeatureSelection(
    spe,
    assay = "counts",
    fam = "binomial",
    sorted = FALSE,
    batch = as.factor(spe$sample_id)
)

pdf(file.path(dir_plots, "binomial_deviance.pdf"))
plot(
    sort(rowData(spe)$binomial_deviance, decreasing = T),
    type = "l",
    xlab = "ranked genes",
    ylab = "binomial deviance",
    main = "Feature Selection with Deviance"
)
abline(v = 2000, lty = 2, col = "red")
dev.off()

message(Sys.time(), " - Running nullResiduals()")
spe <- nullResiduals( # default params
    spe,
    assay = "counts",
    fam = "binomial",
    type = "deviance"
)


hdgs.hb.2000 <-
    rownames(spe)[order(rowData(spe)$binomial_deviance, decreasing = TRUE)][1:2000]
hdgs.hb.5000 <-
    rownames(spe)[order(rowData(spe)$binomial_deviance, decreasing = TRUE)][1:5000]
hdgs.hb.10000 <-
    rownames(spe)[order(rowData(spe)$binomial_deviance, decreasing = TRUE)][1:10000]

save(hdgs.hb.2000,
    hdgs.hb.5000,
    hdgs.hb.10000,
    file = file.path(dir_rdata, "hdgs.hb_2.Rdata")
)

message(Sys.time(), " - Running GLM-PCA")
spe <- runPCA(
    spe,
    exprs_values = "binomial_deviance_residuals",
    subset_row = hdgs.hb.2000,
    ncomponents = num_red_dims,
    name = "GLMPCA_approx",
    BSPARAM = BiocSingular::IrlbaParam()
)

spe <- runPCA(
    spe,
    exprs_values = "binomial_deviance_residuals",
    subset_row = hdgs.hb.5000,
    ncomponents = num_red_dims,
    name = "GLMPCA_approx_5000",
    BSPARAM = BiocSingular::IrlbaParam()
)

spe <- runPCA(
    spe,
    exprs_values = "binomial_deviance_residuals",
    subset_row = hdgs.hb.10000,
    ncomponents = num_red_dims,
    name = "GLMPCA_approx_10000",
    BSPARAM = BiocSingular::IrlbaParam()
)

################################################################################
#   Obtain preliminary clusters based on default GLM-PCA and PCA settings
################################################################################

spe$leiden20_PCA <- clusterCells(spe,
    use.dimred = "PCA",
    BLUSPARAM = SNNGraphParam(
        k = 20,
        cluster.fun = "leiden"
    )
)

spe$leiden20_GLMPCA <- clusterCells(spe,
    use.dimred = "GLMPCA_approx",
    BLUSPARAM = SNNGraphParam(
        k = 20,
        cluster.fun = "leiden"
    )
)


################################################################################
#   Save the processed SPE object
################################################################################

# message(Sys.time(), " - Saving HDF5-backed filtered spe")
# spe = saveHDF5SummarizedExperiment(
#     spe, dir = filtered_hdf5_dir, replace = TRUE
# )
# spe = realize(spe)

message(Sys.time(), " - Saving ordinary filtered spe")
saveRDS(spe, filtered_ordinary_path)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

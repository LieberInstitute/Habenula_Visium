library("ggplot2")
library("SpotSweeper")
library("spatialLIBD")
library("purrr")
library("here")
library("sessioninfo")
library("ggpubr")

## code adapted from https://github.com/LieberInstitute/LFF_spatial_ERC/blob/ad7546fddf48d4047bbc95ab7c980961ac0a3549/code/02_build_spe/04_SpotSweeper.R

plot_dir <- here("plots", "02_build_spe", "03_SpotSweeper")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "02_build_spe", "03_SpotSweeper")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)


set.seed(20241105)
spe <- readRDS(here("processed-data", "02_build_spe", "spe_raw.rds"))

## Remove spots without counts
# if (any(colSums(counts(spe)) == 0)) {
#   message("Removing spots without counts for spe")
#   spe <- spe[, -which(colSums(counts(spe)) == 0)]
#   dim(spe)
# }

# drop out-of-tissue spots
table(spe$in_tissue)
# FALSE  TRUE 
# 35136 44736 

#select the in tissue spots
spe <- spe[, spe$in_tissue]

colnames(colData(spe))

dim(spe)
#[1] 18515 73016

cat("Spots in tissue:", dim(spe)[2], "\n")
#Spots in tissue: 73016

lobstr::obj_size(spe)
# 5.96 GB

#### Spot sweeper ####
# library size
## =========for low sum_umi outliers=====================================
# message(Sys.time(), "Local Outliers - sum_umi")
# spe <- localOutliers(spe,
#                      metric = "sum_umi",
#                      direction = "lower",
#                      log = TRUE
# )
# table(spe$sum_umi_outliers)
# FALSE  TRUE
# 72941    76


## =========for high sum_umi outliers=====================================
# message(Sys.time(), "Local Outliers - high sum_umi")
# spe <- localOutliers(spe,
#                      metric = "sum_umi",
#                      direction = "higher",
#                      log = TRUE
# )
# table(spe$sum_umi_outliers)
# FALSE  TRUE
# 79816    55
# samples <- unique(colData(spe)$sample_id)
# pdf(here(plot_dir, "high_sum_umi_outlier.pdf"), width = 10, height = 8)
# for (sample in samples) {
#     # Generate the plot for the each sample
#     plot <- plotQC(
#         spe = spe,
#         sample_id = "sample_id",
#         sample = sample,
#         metric = "sum_umi",
#         outliers = "sum_umi_outliers",
#         stroke = 0.75
#     )
#
#   # Print the plot to the PDF
#     print(plot)
# }
# dev.off()

## =========for both high and low sum_umi outliers=====================================
# message(Sys.time(), "Local Outliers - high and low sum_umi")
spe <- localOutliers(spe,
                     metric = "sum_umi",
                     direction = "lower",
                     log = TRUE
)
table(spe$sum_umi_outliers)
# FALSE  TRUE 
# 44659    77 


# unique genes
message(Sys.time(), "Local Outliers - sum_gene")
spe <- localOutliers(spe,
                     metric = "sum_gene",
                     direction = "lower",
                     log = TRUE
)
table(spe$sum_gene_outliers)
# FALSE  TRUE 
# 44671    65 

# mitochondrial percent
message(Sys.time(), "Local Outliers - expr_chrM_ratio")
spe <- localOutliers(spe,
                     metric = "expr_chrM_ratio",
                     direction = "higher",
                     log = FALSE
)
table(spe$expr_chrM_ratio_outliers)
# FALSE  TRUE 
# 44696    40 

# combine all outliers into "local_outliers" column
spe$local_outliers <- as.logical(spe$sum_umi_outliers) |
  as.logical(spe$sum_gene_outliers) |
  as.logical(spe$expr_chrM_ratio_outliers)

message("Local Outliers")
table(spe$local_outliers)
# FALSE  TRUE 
# 44611   125

#message("Local Outliers on edge")
# spe$edge_spot
# NULL
#table(spe$local_outliers, spe$edge_spot)

#### find artifacts using SpotSweeper ####
## Only works one sample at a time

unique(spe$sample_id)
# [1] "V13B23-285_A1" "V13B23-285_B1" "V13B23-285_C1" "V13B23-285_D1"
# [5] "V13B23-281_A1" "V13B23-281_B1" "V13B23-281_C1" "V13B23-281_D1"
# [9] "V14F07-340_A1" "V14F07-340_B1" "V14F07-340_C1" "V14F07-340_D1"
# [13] "V13B23-280_A1" "V13B23-280_B1" "V13B23-280_C1" "V13B23-280_D1"
unique(spe$key)
colnames(colData(spe))

message(Sys.time(), " - findArtifact")
artifact_df <- purrr::map_dfr(unique(spe$sample_id), function(samp){
  message(Sys.time(), " - ", samp)
  spe_temp <- findArtifacts(spe[,spe$sample_id == samp],
                            mito_percent = "expr_chrM_ratio",
                            mito_sum = "expr_chrM",
                            n_rings = 5,
                            name = "artifact"
  )

  return(as.data.frame(colData(spe_temp)[,c("sample_id", "key", "artifact")]))
})

head(artifact_df)

## Add artifact to spe
identical(spe$key, artifact_df$key) #TRUE
spe$artifact <- artifact_df$artifact


############## testing one by one sample 

# spe_temp <- findArtifacts(spe[,spe$sample_id == "V13B23-280_C1"],
#                           mito_percent = "expr_chrM_ratio",
#                           mito_sum = "expr_chrM",
#                           n_rings = 5,
#                           name = "artifact"
# )
# # check that "artifact" is now in colData
# colnames(colData(spe_temp))
# # Error: slides  V13B23-280_A1 and V13B23-280_C1 
# # Error in svd(x, nu = 0, nv = k) : infinite or missing values in 'x'

# ## testing
# tmp_lst_samp <- c("V13B23-285_A1", "V13B23-285_B1", "V13B23-285_C1", "V13B23-285_D1", "V13B23-281_A1", "V13B23-281_B1", "V13B23-281_C1", "V13B23-281_D1", "V14F07-340_A1", "V14F07-340_B1", "V14F07-340_C1", "V14F07-340_D1") #"V13B23-280_B1" "V13B23-280_D1"
# length(tmp_lst_samp)
# spe_x <- spe[,spe$sample_id == tmp_lst_samp]
# unique(spe_x$sample_id)

# ## testing removing samples that report issues. All 4 samples for V13B23-280 
#
# message(Sys.time(), " - findArtifact")
# artifact_df <- purrr::map_dfr(unique(spe_x$sample_id), function(samp){
#   message(Sys.time(), " - ", samp)
#   spe_temp <- findArtifacts(spe_x[,spe_x$sample_id == samp],
#                             mito_percent = "expr_chrM_ratio",
#                             mito_sum = "expr_chrM",
#                             n_rings = 5,
#                             name = "artifact"
#   )
# 
#   return(as.data.frame(colData(spe_temp)[,c("sample_id", "key", "artifact")]))
# })

############## END testing 


#### save spot sweeper data ####

spotsweeper_data <- as.data.frame(colData(spe)[,c("sample_id", "key", "array_row", "array_col", "sum_umi_outliers", "sum_gene_outliers", "expr_chrM_ratio_outliers", "local_outliers", "artifact")])
head(spotsweeper_data)
write.csv(spotsweeper_data, file = here(data_dir, "SpotSweeper_data_all_both_dir.csv"))

point_size = 1.1

#### plotting ####

plot_all_spot_sweep <- function(spe, sample = unique(spe$sample_id)[1]){
  
  spe <- spe[,spe$sample_id == sample]
  
  # library size
  p1 <- plotQC(spe, metric = "sum_umi_log", outliers = "sum_umi_outliers", point_size = point_size) +
    ggtitle(paste(sample, "Sum UMI"))
  
  # unique genes
  p2 <- plotQC(spe, metric = "sum_gene_log", outliers = "sum_gene_outliers", point_size = point_size) +
    ggtitle("Sum Genes")
  
  # mitochondrial percent
  p3 <- plotQC(spe, metric = "expr_chrM_ratio", outliers = "expr_chrM_ratio_outliers", point_size = point_size) +
    ggtitle("ChrM Ratio")
  
  # all local outliers
  p4 <- plotQC(spe, metric = "sum_umi_log", outliers = "local_outliers", point_size = point_size, stroke = 0.75) +
    ggtitle("All Local Outliers")
  
  # ## artifact
  # p5 <- plotQC(spe, metric = "sum_umi_log", outliers = "artifact", point_size = point_size, stroke = 0.75) +
  #   ggtitle("Artifact")
  
  # plot_list <- list(p1, p2, p3, p4, p5)
  plot_list <- list(p1, p2, p3, p4)
  ggarrange(
    plotlist = plot_list,
    ncol = 3, nrow = 2,
    common.legend = FALSE)
}

message(Sys.time(), " - Plot SpotSweeper QC for all samples")
pdf(here(plot_dir, "SpotSweeper_ALL_QC_sample.pdf"), width = 12, height = 10)
purrr::map(sort(unique(spe$sample_id)), ~plot_all_spot_sweep(spe = spe, sample = .x))
dev.off()



## Only plot local outliers

plot_all_local_ouliers <- function(spe, sample = unique(spe$sample_id)[1]){
  spe <- spe[,spe$sample_id == sample]
  plt1 <- plotQC(spe, metric = "sum_umi_log", outliers = "local_outliers", 
                 point_size = point_size, stroke = 0.75) +  ggtitle(sample)
}


message("Plot local outliers for all samples")

local_ouliers_plts <- purrr::map(sort(unique(spe$sample_id)), ~ 
                               plot_all_local_ouliers(spe = spe, sample = .x))

length(local_ouliers_plts)
# 16
# local_ouliers_plts[[1]]
pdf(file = here(plot_dir, "SpotSweeper_local_ouliers_v2.pdf"), width = 12, height = 10)

tite_main <- "local_outliers - sum_umi"
plt_main <- ggarrange(
    plotlist = local_ouliers_plts,
    ncol = 4, nrow = 4,
    common.legend = FALSE,
    font.label=list(color="black",size=8)
    )
print(plt_main)

dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


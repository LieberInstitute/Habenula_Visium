
library(SpatialExperiment)
spe <- readRDS("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/ficture_harmony/spe_raw.rds")
head(rowData(spe), 10)
head(colData(spe), 10)
colData(spe)$barcode<-rownames(colData(spe))

library(data.table)
library(dplyr)

repo_dir="/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium"
out_dir=paste0(repo_dir,"/processed-data/10_HD_bin_level/ficture/outputs/all_samples/analysis/nF12.d_12/transcripts_ficture_joined_moved_with_barcodes.tsv.gz")
# read transcripts_ficture_joined
transcripts <- fread(paste0(out_dir), sep = "\t", header = TRUE)

# Compute the number of nonzero factor_K1 values for each barcode
barcode_factor_counts <- rowSums(table(transcripts$barcode, transcripts$factor_K1) > 0)

# Select barcodes that appear in two or more factor_K1 categories
barcodes_with_multiple_factors <- names(barcode_factor_counts[barcode_factor_counts >= 2])

# Display the first 10 barcodes that have multiple factor_K1 assignments
print(head(barcodes_with_multiple_factors, 10)) #character(0) --> Each barcode corresponds to one factor!

dim(table(transcripts$barcode,transcripts$sample_id))

#2555180 

agg_factors <- transcripts[, .(
  factor_K1 = mean(factor_K1, na.rm = TRUE),
  factor_K2 = mean(factor_K2, na.rm = TRUE),
  factor_K3 = mean(factor_K3, na.rm = TRUE),
  factor_P1 = mean(factor_P1, na.rm = TRUE),
  factor_P2 = mean(factor_P2, na.rm = TRUE),
  factor_P3 = mean(factor_P3, na.rm = TRUE)
), by = .(barcode,sample_id)] 

col_data_df <- as.data.frame(colData(spe))
col_data_df <- inner_join(col_data_df, agg_factors, by = c("barcode","sample_id"))
dim(col_data_df)

col_data_df$key<-paste0(col_data_df$barcode,"_",col_data_df$sample_id)
out<-col_data_df[,c("key","factor_K1","factor_K2","factor_K3","factor_P1","factor_P2","factor_P3")]
out<-as.data.frame(out)
print(table(out$factor_K1))
fwrite(out,file="/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/spe_raw/ficture/spe_cluster_all_sample.csv")

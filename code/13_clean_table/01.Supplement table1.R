# Mean ratio markers for cellular, extracellular, multiome (RNA) at each resolution: broad, mid, fine
# Read in all the data

library(data.table)

# cellular:
a<-"/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/09_HD_cell_level/no_secondary/MAGMA/gene_sets"
files <- list.files(a, full.names = TRUE)
files <- files[!grepl("^enrichment_", basename(files))]

dat_list <- lapply(files, fread)
names(dat_list) <- basename(files)
merged_dat <- rbindlist(dat_list, use.names = TRUE, fill = TRUE, idcol = "file")

# extracellular:
a<-"/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/09_HD_cell_level/no_secondary/MAGMA/extracellular/gene_sets"
files <- list.files(a, full.names = TRUE)
files <- files[!grepl("^enrichment_", basename(files))]

dat_list <- lapply(files, fread)
names(dat_list) <- basename(files)
merged_dat_2 <- rbindlist(dat_list, use.names = TRUE, fill = TRUE, idcol = "file")

# multiome:
a<-"/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/10_MAGMA/RNA/gene_sets"
files <- list.files(a, full.names = TRUE)
files <- files[!grepl("^enrichment_", basename(files))]

dat_list <- lapply(files, fread)
names(dat_list) <- basename(files)
merged_dat_3 <- rbindlist(dat_list, use.names = TRUE, fill = TRUE, idcol = "file")


#   Create a final data frame with the relevant columns and some additional metadata
library(data.table)

format_gene_sets <- function(dt, dataset_name) {
  mr_col <- intersect(c("MeanRatio", "mean_ratio"), names(dt))[1]
  if (is.na(mr_col)) stop("No MeanRatio/mean_ratio column found.")
  
  if (!"gene_name" %in% names(dt)) {
    dt[, gene_name := NA_character_]
  }
  
  dt[, .(
    cell_type = set_id,
    cell_type_resolution = sub("\\.tsv$", "", file),
    dataset = dataset_name,
    gene_id = gene_id,
    gene_name = gene_name,
    mean_ratio = get(mr_col)
  )]
}

final_dat <- rbindlist(
  list(
    format_gene_sets(copy(merged_dat),   "cellular"),
    format_gene_sets(copy(merged_dat_2), "extracellular"),
    format_gene_sets(copy(merged_dat_3), "multiome")
  ),
  use.names = TRUE,
  fill = TRUE
)

head(final_dat)

#   Save the final data frame as an RDS file
fwrite(final_dat, "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/13_clean_table/all_gene_sets_merged.csv", sep = "\t")

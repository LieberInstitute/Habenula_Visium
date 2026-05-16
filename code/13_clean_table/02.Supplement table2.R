# Subset of those markers to the genes contributing significant risk for different GWAS traits
# Read in all the data

library(data.table)

# cellular:
a<-"/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/09_HD_cell_level/no_secondary/MAGMA/top_genes.csv"
cellular<-fread(a)

# extracellular:
a<-"/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/09_HD_cell_level/no_secondary/MAGMA/extracellular/top_genes.csv"
extracellular<-fread(a) 

# multiome:
a <- "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/10_MAGMA/RNA/top_genes.csv"
multiome <- fread(a)

#  Save the data

combined <- data.table::rbindlist(
  list(
    cellular = cellular,
    extracellular = extracellular,
    multiome = multiome
  ),
  use.names = TRUE,
  fill = TRUE,
  idcol = "dataset"
)

write.csv(
  combined,
  "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/13_clean_table/supplementary_table2.csv",
  row.names = FALSE
)
## Manisha Barse, Jan 2026

# script to get top enrichment markers of cluster 4 in "09_HD_cell_level/new_samples2/registration_banksy/modeling_results/lambda0_2/1_7.rds" - arranged by FDR" 
library("spatialLIBD")
library("tidyverse")
library("here")

data_dir <- here('processed-data', '09_HD_cell_level', 'new_samples2', 'registration_banksy')
#if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

modeling_results_1_7 <- readRDS(here(data_dir, 'modeling_results', 'lambda0_2', '1_7.rds'))
sig_genes_1_7 <- sig_genes_extract(
    n =20, ## number of top enrichment genes, default=10
    modeling_results = modeling_results_1_7,
    model_type = names(modeling_results_1_7)[2], #enrichment marker genes
    reverse = FALSE,
    sce_layer = fetch_data(type = "sce_layer")
)

head(sig_genes_1_7)
#   top model_type test      gene     stat         pval         fdr gene_index
# 1   1 enrichment   X1 LINC02447 5.383964 2.898197e-07 0.003809551       5097
# 2   2 enrichment   X1   RAPGEF6 5.264336 5.020494e-07 0.003809551       6375
# 3   3 enrichment   X1   ZMYND12 5.101405 1.048076e-06 0.005301868        642
# 4   4 enrichment   X1   SOX1-OT 4.835093 3.381238e-06 0.012828416      15067
# 5   5 enrichment   X1    SORCS2 4.384841 2.228248e-05 0.067631783       5098
# 6   6 enrichment   X1     STAT1 4.150963 5.646056e-05 0.129916378       3389
#      logFC         ensembl
# 1 4.177787 ENSG00000245468
# 2 2.454305 ENSG00000158987
# 3 2.674914 ENSG00000066185
# 4 3.550653 ENSG00000224243
# 5 3.084440 ENSG00000184985
# 6 3.913155 ENSG00000115415

unique(sig_genes_1_7$test)

### get top 20 enrichment genes for cluster4, arranged by FDR
genes_cluster4 <- sig_genes_1_7 |> filter(test == "X4") |> arrange(fdr) #cluster4

## function to save the sig genes
save_gene_names_to_csv <- function(data, filename) {
  file_path <- file.path(data_dir, filename)
  write.csv(data, file = file_path, row.names = FALSE)
}
save_gene_names_to_csv(genes_cluster4, "top20_enrichment_marker_genes_1_7_cluster4.csv")


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
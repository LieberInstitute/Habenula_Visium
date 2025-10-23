# Run GO analysis on the top ligand–receptor pairs

# =====================================================================
# Overall top ligand–receptor pairs
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(enrichplot)
library(tidyr)
library(purrr)
library(tibble)
library(readr)
library(ggplot2)

overall_top_pairs<-read.csv("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/LIANA/table/overall_mean_morans_across_donors.csv")
sce <- readLines("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/LIANA/table/universe_genes.txt")


overall_long <- rbind(
  overall_top_pairs[, c("ligand", "mean", "morans")] |> 
    dplyr::rename(gene = ligand),
  overall_top_pairs[, c("receptor", "mean", "morans")] |> 
    dplyr::rename(gene = receptor)
)

top_genes <- overall_long |>
  dplyr::arrange(desc(mean)) |>
  dplyr::slice(1:30) |>
  dplyr::pull(gene) |> unique()

top_morans_genes <- overall_long |>
  dplyr::arrange(desc(morans)) |>
  dplyr::slice(1:30) |>
  dplyr::pull(gene)|> unique()

fdr_cutoff <- 0.05
plot_dir<-"/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/LIANA/GO"
dir.create(plot_dir, showWarnings = FALSE)

for (ont_type in c("BP", "MF", "CC")) {
  go_obj <- enrichGO(
    gene          = top_genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = ont_type,
    universe = sce,
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    readable      = TRUE
  )
    go_res <- as.data.frame(go_obj) %>%
    filter(p.adjust < fdr_cutoff)

  if (nrow(go_res) > 0) {
    pdf(
      file = file.path(
        plot_dir,
        sprintf("overall_means_GO_%s_topgenes.pdf", ont_type)
      ),
      height = 10,
      width = 8
    )
    print(dotplot(go_obj, showCategory = nrow(go_res)))
    dev.off()
  } else {
    message(sprintf("[skip] %s: no significant terms", ont_type))
  }
}

for (ont_type in c("BP", "MF", "CC")) {
  go_obj <- enrichGO(
    gene          = top_morans_genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = ont_type,
    universe = sce,
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    readable      = TRUE
  )
    go_res <- as.data.frame(go_obj) %>%
    filter(p.adjust < fdr_cutoff)

  if (nrow(go_res) > 0) {
    pdf(
      file = file.path(
        plot_dir,
        sprintf("overall_morans_GO_%s_topgenes.pdf", ont_type)
      ),
      height = 10,
      width = 8
    )
    print(dotplot(go_obj, showCategory = nrow(go_res)))
    dev.off()
  } else {
    message(sprintf("[skip] %s: no significant terms", ont_type))
  }
}

# ======================================================================
# cell-type specific top ligand–receptor pairs

cell_top_pairs<-read.csv("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/LIANA/table/celltype_specific_interactions/celltype_specific_interactions_all.csv")

top_genes_per_cell <- cell_top_pairs %>%
  separate(interaction, into = c("gene1", "gene2"), sep = "\\^") %>%
  pivot_longer(cols = c(gene1, gene2), names_to = "gene_pos", values_to = "gene") %>%
  select(cell_type, gene, abs_z) %>%
  distinct() %>% 
  group_by(cell_type) %>%
  arrange(desc(abs_z),gene,.by_group = TRUE) %>%
  slice_head(n = 30) %>%
  ungroup()

gene_list <- top_genes_per_cell %>%
  group_by(cell_type) %>%
  summarise(genes = list(unique(gene))) %>%
  deframe()

str(gene_list, max.level = 1)
sapply(gene_list, length)

for (ont_type in c("BP", "MF", "CC")) {
    go_obj = compareCluster(
        gene_list, fun = "enrichGO", universe = rownames(sce),
        OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH",
        pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE, keyType = "SYMBOL"
    )
    
    go_obj@compareClusterResult = go_obj@compareClusterResult |>
        filter(p.adjust < fdr_cutoff, ONTOLOGY == ont_type)

    if(nrow(go_obj@compareClusterResult) > 0) {
        pdf(
            file.path(
                plot_dir,
                sprintf(
                    'GO_cell_type_%s_30.pdf',
                    ont_type
                )
            ),
            height = as.integer(
                round(min(18, 2 + nrow(go_obj@compareClusterResult) * 1.5))),
            width = 10
        )
        print(dotplot(go_obj, showCategory = 10)+ theme(
        axis.text.x = element_text(
        angle = 45,     
        hjust = 1,      
        vjust = 1)))
        dev.off()
    }
}

# ======================================================================
# NMF


NMF_top_pairs<-read.csv("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/LIANA/table/NMF/NMF_H_loadings.csv")

# distinct pairs for each factor
weights <- as.matrix(NMF_top_pairs[, -1])
rownames(weights) <- NMF_top_pairs$index

# weights: rows = pairs, columns = factors
calc_z_scores <- function(weights) {
  z_scores <- matrix(NA, nrow = nrow(weights), ncol = ncol(weights))
  rownames(z_scores) <- rownames(weights)
  colnames(z_scores) <- colnames(weights)
  
  for (i in 1:nrow(weights)) {
    for (j in 1:ncol(weights)) {
      other_factors <- weights[i, -j]
      z_scores[i, j] <- (weights[i, j] - mean(other_factors)) / sd(other_factors)
    }
  }
  return(z_scores)
}

# Example usage:
z_scores <- as.data.frame(calc_z_scores(weights))
z_scores$index <- rownames(z_scores)

nmf_long <- z_scores %>%
  separate(index, into = c("gene1", "gene2"), sep = "\\^") %>%
  pivot_longer(cols = starts_with("Factor"), names_to = "factor", values_to = "loading") %>%
  pivot_longer(cols = c(gene1, gene2), names_to = "role", values_to = "gene") %>%
  select(-role)

top_genes_per_factor <- nmf_long %>%
  group_by(factor) %>%
  arrange(desc(loading), .by_group = TRUE) %>%
  slice_head(n = 30) %>%
  ungroup()

gene_list <- top_genes_per_factor %>%
  group_by(factor) %>%
  summarise(genes = list(unique(gene))) %>%
  deframe()

str(gene_list, max.level = 1)
sapply(gene_list, length)
head(gene_list$Factor1)

for (ont_type in c("BP", "MF", "CC")) {
    go_obj = compareCluster(
        gene_list, fun = "enrichGO", universe = rownames(sce),
        OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH",
        pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE, keyType = "SYMBOL"
    )

    go_obj@compareClusterResult = go_obj@compareClusterResult |>
        filter(p.adjust < fdr_cutoff, ONTOLOGY == ont_type)

    if(nrow(go_obj@compareClusterResult) > 0) {
        pdf(
            file.path(
                plot_dir,
                sprintf(
                    'GO_NMF_%s_30.pdf',
                    ont_type
                )
            ),
            height = as.integer(
                round(min(15, 2 + nrow(go_obj@compareClusterResult) * 0.7))),
            width = 10
        )
        print(dotplot(go_obj, showCategory = 10)+ theme(
        axis.text.x = element_text(
        angle = 45,     
        hjust = 1,      
        vjust = 1)))
        dev.off()
    }
}

# ======================================================================
# NMF + cell-types

cell_types<-unique(cell_top_pairs$cell_type)
cell_types[4]<-"Endo_Microglia"

for (cell_type in cell_types){
        data<-read.csv(paste0("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/LIANA/table/NMF/NMF_H_loadings_",cell_type,".csv"))
        
        # distinct pairs for each factor
        weights <- as.matrix(data[, -1])
        rownames(weights) <- data$index

        # weights: rows = pairs, columns = factors
        calc_z_scores <- function(weights) {
          z_scores <- matrix(NA, nrow = nrow(weights), ncol = ncol(weights))
          rownames(z_scores) <- rownames(weights)
          colnames(z_scores) <- colnames(weights)
          
          for (i in 1:nrow(weights)) {
            for (j in 1:ncol(weights)) {
              other_factors <- weights[i, -j]
              z_scores[i, j] <- (weights[i, j] - mean(other_factors)) / sd(other_factors)
            }
          }
          return(z_scores)
        }

        # Example usage:
        z_scores <- as.data.frame(calc_z_scores(weights))
        z_scores$index <- rownames(z_scores)
               
        nmf_long <- z_scores %>%
          separate(index, into = c("gene1", "gene2"), sep = "\\^") %>%
          pivot_longer(cols = starts_with("Factor"), names_to = "factor", values_to = "loading") %>%
          pivot_longer(cols = c(gene1, gene2), names_to = "role", values_to = "gene") %>%
          select(-role)
            
        top_genes_per_factor <- nmf_long %>%
        group_by(factor) %>%
        arrange(desc(loading), .by_group = TRUE) %>%
        slice_head(n = 30) %>%
        ungroup()

        gene_list <- top_genes_per_factor %>%
        group_by(factor) %>%
        summarise(genes = list(unique(gene))) %>%
        deframe()

for (ont_type in c("BP", "MF", "CC")) {
    go_obj = compareCluster(
        gene_list, fun = "enrichGO", universe = rownames(sce),
        OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH",
        pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE, keyType = "SYMBOL"
    )

    go_obj@compareClusterResult = go_obj@compareClusterResult |>
        filter(p.adjust < fdr_cutoff, ONTOLOGY == ont_type)
    
    if(nrow(go_obj@compareClusterResult) > 0) {
        pdf(
            file.path(
                plot_dir,
                sprintf(
                    'GO_NMF_%s_%s.pdf',
                    cell_type, ont_type 
                )
            ),
            height = as.integer(
                round(min(15, 2 + nrow(go_obj@compareClusterResult) * 0.7))),
            width = 10
        )
        print(dotplot(go_obj, showCategory = 10)+ theme(
        axis.text.x = element_text(
        angle = 45,     
        hjust = 1,      
        vjust = 1)))
        dev.off()
    }
}}

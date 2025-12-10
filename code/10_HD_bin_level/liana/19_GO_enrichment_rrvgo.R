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
library(rrvgo)
library(GOSemSim)
semData_BP <- godata('org.Hs.eg.db', ont = "BP")
semData_CC <- godata('org.Hs.eg.db', ont = "CC")
semData_MF <- godata('org.Hs.eg.db', ont = "MF")

fdr_cutoff <- 0.05
num_go_terms <- 20
task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

if(task_id==1){
plot_dir<-"/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/new_samples/liana/GO/cellular"
dir.create(plot_dir, showWarnings = FALSE)
overall_top_pairs<-read.csv("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/new_samples/liana/table/overall_mean_morans_across_donors.csv")
cell_top_pairs<-read.csv("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/new_samples/liana/table/celltype_specific_interactions/celltype_specific_interactions_all.csv")
NMF_top_pairs<-read.csv("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/new_samples/liana/table/NMF/NMF_H_loadings.csv")
sce <- readLines("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/new_samples/liana/table/NMF/universe_genes.txt")

} else{
plot_dir<-"/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/new_samples/liana/GO/extracellular"
dir.create(plot_dir, showWarnings = FALSE)
overall_top_pairs<-read.csv("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/new_samples/liana/table/overall_mean_morans_across_donors_extracellular.csv")
cell_top_pairs<-read.csv("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/new_samples/liana/table/celltype_specific_interactions_extracellular/celltype_specific_interactions_all.csv")
NMF_top_pairs<-read.csv("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/new_samples/liana/table/NMF_extracellular/NMF_H_loadings.csv")
sce <- readLines("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/new_samples/liana/table/NMF_extracellular/universe_genes.txt")
}


overall_long <- rbind(
  overall_top_pairs[, c("ligand", "mean", "morans")] |> 
    dplyr::rename(gene = ligand),
  overall_top_pairs[, c("receptor", "mean", "morans")] |> 
    dplyr::rename(gene = receptor)
)

top_genes <- overall_long |>
  dplyr::arrange(desc(mean)) |>
  dplyr::slice(1:50) |>
  dplyr::pull(gene) |> unique()

top_morans_genes <- overall_long |>
  dplyr::arrange(desc(morans)) |>
  dplyr::slice(1:50) |>
  dplyr::pull(gene)|> unique()

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
    filter(p.adjust < fdr_cutoff) %>%
    arrange(p.adjust) %>%
    slice_head(n = num_go_terms) %>%
    select(ID, Description, p.adjust)

  if (nrow(go_res) > 0) {
    scores <- setNames(-log10(go_res$p.adjust), go_res$ID)
    simMatrix <- calculateSimMatrix(
      go_res$ID,
      orgdb = "org.Hs.eg.db",
      ont   = ont_type,
      method = "Rel",
      semdata = if (ont_type == "BP") {
        semData_BP
      } else if (ont_type == "MF") {
        semData_MF
      } else if (ont_type == "CC") {
        semData_CC
      } else {
        NULL
      }
    )
    reducedTerms <- reduceSimMatrix(
    simMatrix,
    scores = scores,     
    threshold = 0.7,     
    orgdb = "org.Hs.eg.db")

    pdf(
      file = file.path(
        plot_dir,
        sprintf("overall_means_GO_%s_treemap.pdf", ont_type)
      ),
      height = 10,
      width = 8
    )
    invisible(treemapPlot(reducedTerms))
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
    filter(p.adjust < fdr_cutoff) %>%
    arrange(p.adjust) %>%
    slice_head(n = num_go_terms) %>%
    select(ID, Description, p.adjust)

  if (nrow(go_res) > 0) {
    scores <- setNames(-log10(go_res$p.adjust), go_res$ID)
    simMatrix <- calculateSimMatrix(
      go_res$ID,
      orgdb = "org.Hs.eg.db",
      ont   = ont_type,
      method = "Rel",
      semdata = if (ont_type == "BP") {
        semData_BP
      } else if (ont_type == "MF") {
        semData_MF
      } else if (ont_type == "CC") {
        semData_CC
      } else {
        NULL
      }
    )
    reducedTerms <- reduceSimMatrix(
    simMatrix,
    scores = scores,     
    threshold = 0.7,     
    orgdb = "org.Hs.eg.db")

    pdf(
      file = file.path(
        plot_dir,
        sprintf("overall_morans_GO_%s_treemap.pdf", ont_type)
      ),
      height = 10,
      width = 8
    )
    invisible(treemapPlot(reducedTerms))
    dev.off()
  } else {
    message(sprintf("[skip] %s: no significant terms", ont_type))
  }
}

# ======================================================================
# cell-type specific top ligand–receptor pairs

cell_top_pairs$cell_type <- gsub("/", "_", cell_top_pairs$cell_type)

top_genes_per_cell <- cell_top_pairs %>%
  separate(interaction, into = c("gene1", "gene2"), sep = "\\^") %>%
  pivot_longer(cols = c(gene1, gene2), names_to = "gene_pos", values_to = "gene") %>%
  select(cell_type, gene, z_spec) %>%
  distinct() %>% 
  group_by(cell_type) %>%
  arrange(desc(z_spec),gene,.by_group = TRUE) %>%
  filter(z_spec>=2) %>%
  ungroup()

gene_list <- top_genes_per_cell %>%
  group_by(cell_type) %>%
  summarise(genes = list(unique(gene))) %>%
  deframe()

str(gene_list, max.level = 1)
sapply(gene_list, length)

go_obj = compareCluster(
        gene_list, fun = "enrichGO", universe = rownames(sce),
        OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH",
        pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE, keyType = "SYMBOL"
    )

for (ont_type in c("BP", "MF", "CC")) {   
  for (cluster in names(gene_list)) { 
    go_res <- as.data.frame(go_obj) %>%
    filter(p.adjust < fdr_cutoff, ONTOLOGY == ont_type, Cluster == cluster) %>%
    arrange(p.adjust) %>%
    slice_head(n = num_go_terms) %>%
    select(ID, Description, p.adjust)

    if (nrow(go_res) > 0) {
    scores <- setNames(-log10(go_res$p.adjust), go_res$ID)
    simMatrix <- calculateSimMatrix(
      go_res$ID,
      orgdb = "org.Hs.eg.db",
      ont   = ont_type,
      method = "Rel",
      semdata = if (ont_type == "BP") {
        semData_BP
      } else if (ont_type == "MF") {
        semData_MF
      } else if (ont_type == "CC") {
        semData_CC
      } else {
        NULL
      }
    )
    reducedTerms <- reduceSimMatrix(
    simMatrix,
    scores = scores,     
    threshold = 0.7,     
    orgdb = "org.Hs.eg.db")

    pdf(
      file = file.path(
        plot_dir,
        sprintf("GO_cell_type_%s_z2_%s_treemap.pdf", ont_type, cluster)
      ),
      height = 20,
      width = 10
    )
    invisible(treemapPlot(reducedTerms))
    dev.off()
  } else {
    message(sprintf("[skip] %s: no significant terms", ont_type))
  }
}}

# ======================================================================
# NMF

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

# z-score >=2

top_genes_per_factor <- nmf_long %>%
  group_by(factor) %>%
  arrange(desc(loading), .by_group = TRUE) %>%
  filter(loading>=2) %>%
  ungroup()

gene_list <- top_genes_per_factor %>%
  group_by(factor) %>%
  summarise(genes = list(unique(gene))) %>%
  deframe()

str(gene_list, max.level = 1)
sapply(gene_list, length)
head(gene_list$Factor1)

go_obj = compareCluster(
        gene_list, fun = "enrichGO", universe = rownames(sce),
        OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH",
        pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE, keyType = "SYMBOL"
    )

for (ont_type in c("BP", "MF", "CC")) {
  for (cluster in names(gene_list)) {
    message(sprintf("Cluster: %s, Ontology: %s", cluster, ont_type))
  go_res <- as.data.frame(go_obj) %>%
    filter(p.adjust < fdr_cutoff, ONTOLOGY == ont_type, Cluster == cluster) %>%
    arrange(p.adjust) %>%
    slice_head(n = num_go_terms) %>%
    select(ID, Description, p.adjust)

    if (nrow(go_res) > 0) {
    scores <- setNames(-log10(go_res$p.adjust), go_res$ID)
    simMatrix <- calculateSimMatrix(
      go_res$ID,
      orgdb = "org.Hs.eg.db",
      ont   = ont_type,
      method = "Rel",
      semdata = if (ont_type == "BP") {
        semData_BP
      } else if (ont_type == "MF") {
        semData_MF
      } else if (ont_type == "CC") {
        semData_CC
      } else {
        NULL
      }
    )

    reducedTerms <- reduceSimMatrix(
    simMatrix,
    scores = scores,     
    threshold = 0.7,     
    orgdb = "org.Hs.eg.db")

    pdf(
          file = file.path(
            plot_dir,
            sprintf("GO_NMF_%s_z2_%s_treemap.pdf", ont_type, cluster)
          ),
          height = 20,
          width = 10
        )
        invisible(treemapPlot(reducedTerms))
        dev.off()
      } else {
        message(sprintf("[skip] %s: no significant terms", ont_type))
      }
    }}


# ======================================================================
# NMF + cell-types

cell_types<-unique(cell_top_pairs$cell_type)

if (task_id==1){
for (cell_type in cell_types){
        data<-read.csv(paste0("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/new_samples/liana/table/NMF/NMF_H_loadings_",cell_type,".csv"))
        
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
        filter(loading>=2) %>%
        ungroup()

        gene_list <- top_genes_per_factor %>%
        group_by(factor) %>%
        summarise(genes = list(unique(gene))) %>%
        deframe()

    go_obj = compareCluster(
        gene_list, fun = "enrichGO", universe = rownames(sce),
        OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH",
        pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE, keyType = "SYMBOL"
    )

for (ont_type in c("BP", "MF", "CC")) {
 for (cluster in names(gene_list)) {
    message(sprintf("Cell type: %s, Factor: %s, Ontology: %s", cell_type, cluster, ont_type))
  
  go_res <- as.data.frame(go_obj) %>%
    filter(p.adjust < fdr_cutoff, ONTOLOGY == ont_type, Cluster == cluster) %>%
    arrange(p.adjust) %>%
    slice_head(n = num_go_terms) %>%
    select(ID, Description, p.adjust)

    if (nrow(go_res) > 0) {
    scores <- setNames(-log10(go_res$p.adjust), go_res$ID)
    simMatrix <- calculateSimMatrix(
      go_res$ID,
      orgdb = "org.Hs.eg.db",
      ont   = ont_type,
      method = "Rel",
      semdata = if (ont_type == "BP") {
        semData_BP
      } else if (ont_type == "MF") {
        semData_MF
      } else if (ont_type == "CC") {
        semData_CC
      } else {
        NULL
      }
    )

    reducedTerms <- reduceSimMatrix(
    simMatrix,
    scores = scores,     
    threshold = 0.7,     
    orgdb = "org.Hs.eg.db")


    pdf(
          file = file.path(
            plot_dir,
            sprintf("GO_NMF_%s_%s_z2_%s_treemap.pdf", cell_type, ont_type, cluster)
          ),
          height = 20,
          width = 10
        )
        invisible(treemapPlot(reducedTerms))
        dev.off()
      } else {
        message(sprintf("[skip] %s: no significant terms", ont_type))
      }
    }
}}} else{
  for (cell_type in cell_types){
        data<-read.csv(paste0("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/new_samples/liana/table/NMF_extracellular/NMF_H_loadings_",cell_type,".csv"))
        
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
        filter(loading>=2) %>%
        ungroup()

        gene_list <- top_genes_per_factor %>%
        group_by(factor) %>%
        summarise(genes = list(unique(gene))) %>%
        deframe()

    go_obj = compareCluster(
        gene_list, fun = "enrichGO", universe = rownames(sce),
        OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH",
        pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE, keyType = "SYMBOL"
    )

for (ont_type in c("BP", "MF", "CC")) {
 for (cluster in names(gene_list)) {
    message(sprintf("Cell type: %s, Factor: %s, Ontology: %s", cell_type, cluster, ont_type))
  
  go_res <- as.data.frame(go_obj) %>%
    filter(p.adjust < fdr_cutoff, ONTOLOGY == ont_type, Cluster == cluster) %>%
    arrange(p.adjust) %>%
    slice_head(n = num_go_terms) %>%
    select(ID, Description, p.adjust)

    if (nrow(go_res) > 0) {
    scores <- setNames(-log10(go_res$p.adjust), go_res$ID)
    simMatrix <- calculateSimMatrix(
      go_res$ID,
      orgdb = "org.Hs.eg.db",
      ont   = ont_type,
      method = "Rel",
      semdata = if (ont_type == "BP") {
        semData_BP
      } else if (ont_type == "MF") {
        semData_MF
      } else if (ont_type == "CC") {
        semData_CC
      } else {
        NULL
      }
    )

    reducedTerms <- reduceSimMatrix(
    simMatrix,
    scores = scores,     
    threshold = 0.7,     
    orgdb = "org.Hs.eg.db")


    pdf(
          file = file.path(
            plot_dir,
            sprintf("GO_NMF_%s_%s_z2_%s_treemap.pdf", cell_type, ont_type, cluster)
          ),
          height = 20,
          width = 10
        )
        invisible(treemapPlot(reducedTerms))
        dev.off()
      } else {
        message(sprintf("[skip] %s: no significant terms", ont_type))
      }
    }
}}
}

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
library(here)
library(sessioninfo)

semData_BP <- godata(annoDb = "org.Hs.eg.db", ont = "BP")
semData_CC <- godata(annoDb = "org.Hs.eg.db", ont = "CC")
semData_MF <- godata(annoDb = "org.Hs.eg.db", ont = "MF")

fdr_cutoff <- 0.05
num_go_terms <- 50
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

if (task_id == 1) {
    plot_dir <- here(
        'processed-data', '10_HD_bin_level', 'no_secondary', 'liana', 'GO', 'cellular'
    )
    dir.create(plot_dir, showWarnings = FALSE)
    overall_top_pairs <- read.csv("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/table/overall_mean_morans_across_donors.csv")
    cell_top_pairs <- read.csv("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/table/celltype_specific_interactions/celltype_specific_interactions_all.csv")
    NMF_top_pairs <- read.csv("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/table/NMF/NMF_H_loadings.csv")
    sce <- readLines("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/table/NMF/universe_genes.txt")
    nmf_cell_loadings_path <- paste0("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/table/NMF/NMF_H_loadings_%s.csv")
} else {
    plot_dir <- "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/GO/extracellular"
    dir.create(plot_dir, showWarnings = FALSE)
    overall_top_pairs <- read.csv("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/table/overall_mean_morans_across_donors_extracellular.csv")
    cell_top_pairs <- read.csv("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/table/celltype_specific_interactions_extracellular/celltype_specific_interactions_all.csv")
    NMF_top_pairs <- read.csv("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/table/NMF_extracellular/NMF_H_loadings.csv")
    sce <- readLines("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/table/NMF_extracellular/universe_genes.txt")
    nmf_cell_loadings_path <- paste0("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/table/NMF_extracellular/NMF_H_loadings_%s.csv")
}

#-------------------------------------------------------------------------------
#   Functions
#-------------------------------------------------------------------------------

GO_wrapper <- function(gene_list, plot_prefix) {

  # ---- helper: GOID -> Term ----
  getGoTerm <- function(go_ids) {
    if (!requireNamespace("GO.db", quietly = TRUE) ||
        !requireNamespace("AnnotationDbi", quietly = TRUE)) {
      return(rep(NA_character_, length(go_ids)))
    }
    tt <- suppressMessages(AnnotationDbi::Term(GO.db::GOTERM[go_ids]))
    out <- as.character(tt)

    # keep order / length aligned with go_ids
    out2 <- rep(NA_character_, length(go_ids))
    names(out2) <- go_ids
    out2[names(tt)] <- out
    unname(out2[go_ids])
  }

  # ---- 0) accept vector OR list(length==1) ----
  if (is.list(gene_list)) {
    if (length(gene_list) != 1) {
      stop("GO_wrapper(): please call this function per cluster (gene vector). gene_list is a list with length > 1.")
    }
    gene_list <- gene_list[[1]]
  }
  genes <- unique(na.omit(as.character(gene_list)))
  if (length(genes) == 0) stop("GO_wrapper(): gene_list is empty after cleaning")

  # universe safety
  universe_genes <- unique(na.omit(as.character(sce)))

  # ---- 1) run enrichGO ONCE (ALL) ----
  ego_all <- suppressMessages(
    enrichGO(
      gene          = genes,
      universe      = universe_genes,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff  = 1,
      qvalueCutoff  = 1,
      readable      = TRUE
    )
  )

  go_df <- as.data.frame(ego_all)
  if (nrow(go_df) == 0) {
    message(sprintf("[skip] %s: enrichGO returned 0 rows", plot_prefix))
    # still write empty files so pipeline won't break
    write.csv(data.frame(), file.path(plot_dir, sprintf("%s_GO_raw_top%d.csv", plot_prefix, num_go_terms)), row.names = FALSE)
    write.csv(data.frame(), file.path(plot_dir, sprintf("%s_GO_rrvgo_parent_child.csv", plot_prefix)), row.names = FALSE)
    return(invisible(ego_all))
  }

  # ---- 2) RAW table: GO of top LR genes by semantic category ----
  raw_tbl <- go_df %>%
    dplyr::filter(ONTOLOGY %in% c("BP","MF","CC"), p.adjust < fdr_cutoff) %>%
    dplyr::group_by(ONTOLOGY) %>%
    dplyr::arrange(p.adjust, .by_group = TRUE) %>%
    dplyr::slice_head(n = num_go_terms) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      cell_type = plot_prefix,            # 你外面按 cluster 调用，所以这里用 plot_prefix 标记
      semantic_category = ONTOLOGY,
      go_id = ID,
      go_term = Description,
      p.adjust,
      Count,
      geneID
    )

  write.csv(
    raw_tbl,
    file = file.path(plot_dir, sprintf("%s_GO_raw_top%d.csv", plot_prefix, num_go_terms)),
    row.names = FALSE
  )

  # ---- 3) rrvgo parent/child table + treemap per ontology ----
  parent_child_all <- list()

  for (ont_type in c("BP", "MF", "CC")) {

    go_res <- go_df %>%
      dplyr::filter(ONTOLOGY == ont_type, p.adjust < fdr_cutoff) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = num_go_terms) %>%
      dplyr::select(ID, Description, p.adjust)

    if (nrow(go_res) < 2) {
      message(sprintf("[skip] %s | %s: <2 significant terms", plot_prefix, ont_type))
      next
    }

    scores <- setNames(-log10(go_res$p.adjust), go_res$ID)

    simMatrix <- calculateSimMatrix(
      go_res$ID,
      orgdb   = "org.Hs.eg.db",
      ont     = ont_type,
      method  = "Rel",
      semdata = switch(
        ont_type,
        BP = semData_BP,
        MF = semData_MF,
        CC = semData_CC
      )
    )

    reducedTerms <- reduceSimMatrix(
      simMatrix,
      scores   = scores,
      threshold = 0.7,
      orgdb    = "org.Hs.eg.db"
    )

    reduced_df <- as.data.frame(reducedTerms)

    parent_col <- dplyr::case_when(
      "parent"  %in% colnames(reduced_df) ~ "parent",
      "reducer" %in% colnames(reduced_df) ~ "reducer",
      TRUE ~ NA_character_
    )

    pc <- reduced_df %>%
      dplyr::mutate(
        cell_type = plot_prefix,
        semantic_category = ont_type,
        child_go   = .data$go,
        child_term = getGoTerm(.data$go),
        parent_go  = if (is.na(parent_col)) .data$go else dplyr::if_else(
          is.na(.data[[parent_col]]), .data$go, as.character(.data[[parent_col]])
        ),
        parent_term = getGoTerm(.data$parent_go)
      ) %>%
      dplyr::select(
        cell_type, semantic_category,
        parent_go, parent_term,
        child_go, child_term,
        dplyr::any_of("score")
      )

    parent_child_all[[ont_type]] <- pc

    # treemap
    pdf(file = file.path(plot_dir, sprintf("%s_%s_treemap.pdf", plot_prefix, ont_type)))
    invisible(treemapPlot(reduced_df))
    dev.off()
  }

  parent_child_tbl <- dplyr::bind_rows(parent_child_all)

  write.csv(
    parent_child_tbl,
    file = file.path(plot_dir, sprintf("%s_GO_rrvgo_parent_child.csv", plot_prefix)),
    row.names = FALSE
  )

  invisible(ego_all)
}


NMF_cell_gene_list <- function(cell_type) {
    data <- read.csv(sprintf(nmf_cell_loadings_path, cell_type))

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
        filter(loading >= 2) %>%
        ungroup()

    gene_list <- top_genes_per_factor %>%
        group_by(factor) %>%
        summarise(genes = list(unique(gene))) %>%
        deframe()
    
    return(gene_list)
}

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

#-------------------------------------------------------------------------------
#   Main
#-------------------------------------------------------------------------------

overall_long <- rbind(
    overall_top_pairs[, c("ligand", "mean", "morans")] |>
        dplyr::rename(gene = ligand),
    overall_top_pairs[, c("receptor", "mean", "morans")] |>
        dplyr::rename(gene = receptor)
)

top_genes <- overall_long |>
    dplyr::arrange(desc(mean)) |>
    dplyr::slice(1:50) |>
    dplyr::pull(gene) |>
    unique()

top_morans_genes <- overall_long |>
    dplyr::arrange(desc(morans)) |>
    dplyr::slice(1:50) |>
    dplyr::pull(gene) |>
    unique()

GO_wrapper(top_genes, 'overall_means_GO')
GO_wrapper(top_morans_genes, 'overall_morans_GO')

# ======================================================================
# cell-type specific top ligand–receptor pairs

cell_top_pairs$cell_type <- gsub("/", "_", cell_top_pairs$cell_type)

top_genes_per_cell <- cell_top_pairs %>%
    separate(interaction, into = c("gene1", "gene2"), sep = "\\^") %>%
    pivot_longer(cols = c(gene1, gene2), names_to = "gene_pos", values_to = "gene") %>%
    select(cell_type, gene, z_spec) %>%
    distinct() %>%
    group_by(cell_type) %>%
    arrange(desc(z_spec), gene, .by_group = TRUE) %>%
    filter(z_spec >= 2) %>%
    ungroup()

gene_list <- top_genes_per_cell %>%
    group_by(cell_type) %>%
    summarise(genes = list(unique(gene))) %>%
    deframe()

str(gene_list, max.level = 1)
sapply(gene_list, length)

for (cluster in names(gene_list)) {
    GO_wrapper(gene_list[[cluster]], sprintf("GO_cell_type_%s_z2", cluster))
}

# ======================================================================
# NMF

# distinct pairs for each factor
weights <- as.matrix(NMF_top_pairs[, -1])
rownames(weights) <- NMF_top_pairs$index

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
    filter(loading >= 2) %>%
    ungroup()

gene_list <- top_genes_per_factor %>%
    group_by(factor) %>%
    summarise(genes = list(unique(gene))) %>%
    deframe()

str(gene_list, max.level = 1)
sapply(gene_list, length)
head(gene_list$Factor1)

for (cluster in names(gene_list)) {
    GO_wrapper(gene_list[[cluster]], sprintf("GO_NMF_%s_z2", cluster))
}

# ======================================================================
# NMF + cell-types

cell_types <- unique(cell_top_pairs$cell_type)

for (cell_type in cell_types) {
    gene_list <- NMF_cell_gene_list(cell_type)

    for (cluster in names(gene_list)) {
        GO_wrapper(
            gene_list[[cluster]], sprintf("GO_NMF_%s_%s_z2", cell_type, cluster)
        )
    }
}

session_info()

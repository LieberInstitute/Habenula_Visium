suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
  library(stringr); library(purrr); library(tibble)
  library(ComplexHeatmap); library(circlize); library(grid)
  library(viridisLite)
})

# ============================== #
# 0) Global settings
# ============================== #
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(task_id) || task_id == 1) {
  dir_in  <- "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/GO/cellular"
} else {
  dir_in  <- "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/GO/extracellular"
}
pattern    <- "^GO_cell_type_.*_GO_rrvgo_parent_child\\.csv$"
out_prefix <- file.path(dir_in, "GO_celltype_parentChild")

# Tunables
min_cells            <- 1           # child term must appear (>0) in >= N cells
top_per_parent       <- 8           # keep top-N child terms per parent
top_global_cap       <- 300         # global cap on child terms
value_transform_mode <- "none"      # "none" or "row_z"
wrap_width           <- 28          # wrap width for child-term labels
row_name_max_cm      <- 12          # max width for row names
cap_quantile         <- 0.99        # color cap quantile
cluster_cols         <- TRUE        # cluster columns to get order
show_col_dend        <- FALSE       # hide dendrogram (keep order)
family_title_fs      <- 7           # family (slice) title font size on the right

# Margins / spacing
pad_top_mm    <- 2
pad_right_mm  <- 2
pad_bottom_mm <- 2
pad_left_mm   <- 2                  # outer left page margin (kept small)
left_spacer_mm <- 16                # ✅ inner left white spacer beside the heatmap

# ============================== #
# 1) Helpers
# ============================== #
parse_fname <- function(f){
  b <- basename(f)
  m <- str_match(b, "^GO_cell_type_(.+?)_(z\\d+)_GO_rrvgo_parent_child\\.csv$")
  tibble(
    cell_type_name = ifelse(!is.na(m[,2]), m[,2], NA_character_),
    z_label        = ifelse(!is.na(m[,3]), m[,3], NA_character_)
  )
}

pick_col <- function(df, candidates, regex = NULL){
  cn <- names(df)
  nm <- intersect(candidates, cn)[1]
  if(!is.na(nm)) return(nm)
  if(!is.null(regex)){
    nm2 <- cn[str_detect(tolower(cn), regex)][1]
    if(!is.na(nm2)) return(nm2)
  }
  NA_character_
}

# Prefer -log10(p); fallback to score/NES/ES; else presence (=1)
compute_value <- function(df){
  p_cand  <- c("p.adjust","padj","p_adj","pvalue","p_value","PValue","adj_p","p")
  sc_cand <- c("score","Score","NES","enrichmentScore","ES")
  pc <- pick_col(df, p_cand)
  if(!is.na(pc)) return(-log10(pmax(as.numeric(df[[pc]]), .Machine$double.eps)))
  sc <- pick_col(df, sc_cand)
  if(!is.na(sc)) return(as.numeric(df[[sc]]))
  rep(1, nrow(df))
}

wrap_labels <- function(x, width = 36){
  str_replace_all(str_wrap(x, width = width), "\\s*-\\s*", "-")
}

# Clean cell names for x-axis
clean_cell_name <- function(x){
  x <- as.character(x)
  x[is.na(x)] <- "Unknown"
  m <- str_match(x, "^GO_cell_type_(.+?)_(z\\d+)$")
  x <- ifelse(!is.na(m[,2]), m[,2], x)         # strip "GO_cell_type_" and "_zN"
  x <- str_replace_all(x, "\\s*\\(z\\d+\\)", "")
  x
}

# ============================== #
# 2) Robust single-file reader
# ============================== #
read_one <- function(f){
  info <- tryCatch(file.info(f)$size, error = function(e) NA_real_)
  if (is.na(info) || info == 0) { message("[skip] empty file: ", basename(f)); return(NULL) }
  raw <- tryCatch(suppressMessages(read_csv(f, show_col_types = FALSE)), error = function(e) NULL)
  if (is.null(raw) || nrow(raw) == 0) { message("[skip] read failed or header-only: ", basename(f)); return(NULL) }

  col_cell   <- pick_col(raw, c("cell_type","CellType","celltype"))
  col_sem    <- pick_col(raw, c("semantic_category","ONTOLOGY","ontology","Category"), "onto|semantic|category")
  col_parent <- pick_col(raw, c("parent_term","ParentTerm","parent_desc","ParentDescription","Representative_Term"),
                         "parent.*(term|desc|name)|representative")
  col_child  <- pick_col(raw, c("child_term","ChildTerm","term","Term","Description","name","label"),
                         "term|desc|name|label")
  if (is.na(col_child)) { message("[skip] child_term not found: ", basename(f)); return(NULL) }

  meta <- parse_fname(f)

  tibble(
    .file       = basename(f),
    cell_type   = if(!is.na(col_cell)) raw[[col_cell]] else paste0("GO_cell_type_", meta$cell_type_name[1], "_", meta$z_label[1]),
    z_label     = meta$z_label[1],
    ontology    = if(!is.na(col_sem))    raw[[col_sem]]    else NA_character_,
    parent_term = if(!is.na(col_parent)) raw[[col_parent]] else "Unassigned",
    child_term  = raw[[col_child]],
    value       = compute_value(raw)
  ) %>%
    filter(!is.na(child_term) & trimws(child_term) != "")
}

# ============================== #
# 3) Read all and normalize
# ============================== #
files <- list.files(dir_in, pattern = pattern, full.names = TRUE)
stopifnot(length(files) > 0)

combined_long <- map(seq_along(files), function(i){
  f <- files[i]
  tryCatch({
    df <- read_one(f)
    if (!is.null(df)) message(sprintf("[ok] %d/%d %s -> %d rows", i, length(files), basename(f), nrow(df)))
    df
  }, error = function(e){
    message(sprintf("[skip] error reading %s | %s", basename(f), e$message))
    NULL
  })
}) %>% compact() %>% bind_rows()

combined_long <- combined_long %>%
  mutate(
    ontology = case_when(
      is.na(ontology) ~ NA_character_,
      str_to_upper(ontology) %in% c("BP","BIOLOGICAL_PROCESS") ~ "BP",
      str_to_upper(ontology) %in% c("MF","MOLECULAR_FUNCTION") ~ "MF",
      str_to_upper(ontology) %in% c("CC","CELLULAR_COMPONENT") ~ "CC",
      TRUE ~ ontology
    ),
    cell_type = clean_cell_name(cell_type)
  )
message("[info] total rows: ", nrow(combined_long))

# ============================== #
# 4) Plot function
#    - one family label per block on the right
#    - column clustering for order (dendrogram hidden)
#    - extra left white space (inner spacer) beside the heatmap body
# ============================== #
plot_one_ontology <- function(dat, ont_label, out_prefix,
                              min_cells = 1, top_per_parent = 8, top_global_cap = 300,
                              value_transform_mode = c("none","row_z"),
                              wrap_width = 28, row_name_max_cm = 12,
                              cap_quantile = 0.99,
                              cluster_cols = TRUE, show_col_dend = FALSE,
                              family_title_fs = 7,
                              pad_top_mm = 2, pad_right_mm = 2, pad_bottom_mm = 2, pad_left_mm = 2,
                              left_spacer_mm = 16){

  value_transform_mode <- match.arg(value_transform_mode)

  dat <- dat %>% filter(ontology == ont_label)
  if (nrow(dat) == 0) { message("[info] no rows for ", ont_label); return(invisible(NULL)) }

  dat <- dat %>%
    filter(!is.na(child_term), trimws(child_term) != "",
           !is.na(parent_term), trimws(parent_term) != "") %>%
    group_by(cell_type, child_term, parent_term) %>%
    summarise(value = max(value, na.rm = TRUE), .groups = "drop")

  # representative parent per child (prefer self; else highest mean)
  parent_choice <- dat %>%
    group_by(child_term, parent_term) %>%
    summarise(mean_val = mean(value, na.rm = TRUE),
              n = n(), self = any(child_term == parent_term), .groups = "drop") %>%
    arrange(desc(self), desc(mean_val), desc(n), parent_term) %>%
    group_by(child_term) %>% slice(1) %>% ungroup() %>%
    select(child_term, parent_term)

  dat <- dat %>% select(-parent_term) %>% left_join(parent_choice, by = "child_term")

  # keep terms that appear (>0) in at least min_cells distinct cells
  keep_terms <- dat %>%
    filter(value > 0) %>%
    distinct(child_term, cell_type) %>%
    count(child_term, name = "n_cells") %>%
    filter(n_cells >= min_cells) %>%
    pull(child_term)
  dat <- dat %>% filter(child_term %in% keep_terms)
  if (nrow(dat) == 0) { message("[info] none after min_cells filter for ", ont_label); return(invisible(NULL)) }

  # thin per-parent
  if (is.finite(top_per_parent)) {
    keep_by_parent <- dat %>%
      group_by(parent_term, child_term) %>%
      summarise(score = mean(value, na.rm = TRUE), .groups = "drop_last") %>%
      slice_max(score, n = top_per_parent, with_ties = FALSE) %>%
      pull(child_term)
    dat <- dat %>% filter(child_term %in% keep_by_parent)
  }

  # global cap
  if (is.finite(top_global_cap)) {
    keep_global <- dat %>%
      group_by(child_term) %>%
      summarise(score = sum(value, na.rm = TRUE), .groups = "drop") %>%
      slice_max(score, n = top_global_cap, with_ties = FALSE) %>%
      pull(child_term)
    dat <- dat %>% filter(child_term %in% keep_global)
  }
  if (nrow(dat) == 0) { message("[info] none after thinning for ", ont_label); return(invisible(NULL)) }

  # wide matrix
  mat_wide <- dat %>% select(child_term, cell_type, value, parent_term) %>%
    pivot_wider(names_from = cell_type, values_from = value, values_fill = 0)
  row_meta <- mat_wide %>% select(child_term, parent_term)
  mat <- mat_wide %>% select(-child_term, -parent_term) %>% as.data.frame()
  rownames(mat) <- row_meta$child_term

  # row order: parent block total, then within-block total
  parent_score <- dat %>% group_by(parent_term) %>%
    summarise(parent_total = sum(value, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(parent_total))
  parent_levels <- parent_score$parent_term
  parent_levels <- c(setdiff(parent_levels, c("Unassigned","NA","")),
                     intersect(parent_levels, c("Unassigned","NA","")))
  child_tot <- dat %>% group_by(child_term) %>% summarise(tot = sum(value, na.rm = TRUE), .groups = "drop")
  row_order <- row_meta %>% left_join(child_tot, by = "child_term") %>%
    arrange(factor(parent_term, levels = parent_levels), desc(tot), child_term) %>%
    pull(child_term)
  mat <- mat[row_order, , drop = FALSE]
  row_parent <- factor(row_meta$parent_term[match(rownames(mat), row_meta$child_term)], levels = parent_levels)

  # child-term labels
  child_wrapped <- wrap_labels(rownames(mat), width = wrap_width)

  # values & colors
  if (value_transform_mode == "row_z") {
    mat_plot <- t(scale(t(as.matrix(mat))))
    mat_plot[!is.finite(mat_plot)] <- 0
    mat_plot <- pmin(pmax(mat_plot, -3), 3)
    col_fun <- colorRamp2(c(-3, 0, 3), c("#2C7BB6", "#F7F7F7", "#D7191C"))
    legend_at <- c(-3, -1.5, 0, 1.5, 3)
    legend_title <- paste0("z-score (", ont_label, ")")
  } else {
    mat_plot <- as.matrix(mat)
    v <- as.numeric(mat_plot); v <- v[is.finite(v)]; if (length(v) == 0) v <- c(0, 1)
    vmax <- stats::quantile(v, probs = cap_quantile, na.rm = TRUE)
    col_fun <- colorRamp2(c(0, vmax*0.33, vmax*0.66, vmax), viridis(4))
    legend_at <- round(c(0, vmax*0.33, vmax*0.66, vmax), 2)
    legend_title <- paste0("-log10(p) (", ont_label, ")")
  }

  # sizes & fonts
  rownames(mat_plot) <- child_wrapped
  fs_row <- ifelse(nrow(mat_plot) > 150, 4, ifelse(nrow(mat_plot) > 80, 5, 6))
  fs_col <- ifelse(ncol(mat_plot) > 40, 7, 9)
  sz_w <- max(8, min(24, 4 + ncol(mat_plot) * 0.25))
  sz_h <- max(8, min(36, 4 + nrow(mat_plot) * 0.20))

  # blank left spacer (inner white strip)
  left_spacer <- rowAnnotation(
    .sp = anno_empty(border = FALSE),
    width = unit(left_spacer_mm, "mm")
  )

  # main heatmap
  ht <- Heatmap(
    mat_plot,
    name = paste0("score_", ont_label),
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = cluster_cols,
    show_column_dend = show_col_dend,      # hide top dendrogram but keep order
    column_dend_height = unit(0, "mm"),
    row_split = row_parent,                # family blocks
    row_title_side = "right",
    row_title_rot  = 0,
    row_title_gp   = gpar(fontsize = family_title_fs, fontface = "bold"),
    row_gap  = unit(2, "mm"),
    show_row_names = TRUE,
    row_names_side = "right",
    row_names_gp   = gpar(fontsize = fs_row),
    row_names_max_width = unit(row_name_max_cm, "cm"),
    show_column_names = TRUE,
    column_names_gp   = gpar(fontsize = fs_col),
    column_names_rot  = 45,
    heatmap_legend_param = list(
      title = legend_title,
      at = legend_at,
      title_gp = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)
    )
  )

  # combine spacer + heatmap
  htlist <- left_spacer + ht

  # # export
  # out_pdf <- paste0(out_prefix, "_heatmap_", ont_label, ".pdf")
  # pdf(out_pdf, width = sz_w, height = sz_h)
  # draw(
  #   htlist,
  #   heatmap_legend_side = "right",
  #   padding = unit.c(                           # outer page margins
  #     unit(pad_top_mm, "mm"),
  #     unit(pad_right_mm, "mm"),
  #     unit(pad_bottom_mm, "mm"),
  #     unit(pad_left_mm, "mm")
  #   )
  # )
  # dev.off()
  # message("[ok] heatmap saved: ", out_pdf)

  # export (PNG)
  out_png <- paste0(out_prefix, "_heatmap_", ont_label, ".png")
  png(
    filename = out_png,
    width = sz_w, height = sz_h, 
    units = "in", res = 300,
    type = "cairo-png", bg = "white"
  )
  draw(
    htlist,
    heatmap_legend_side = "right",
    padding = unit.c(
      unit(pad_top_mm, "mm"),
      unit(pad_right_mm, "mm"),
      unit(pad_bottom_mm, "mm"),
      unit(pad_left_mm, "mm")
    )
  )
  dev.off()
  message("[ok] heatmap saved: ", out_png)


  out_csv <- paste0(out_prefix, "_matrix_", ont_label, ".csv")
  mat_out <- as.data.frame(mat) %>%
    rownames_to_column("GO_term") %>%
    left_join(row_meta, by = c("GO_term" = "child_term")) %>%
    relocate(GO_term, parent_term)
  write_csv(mat_out, out_csv)
  message("[ok] matrix saved: ", out_csv)

  invisible(list(mat = mat, row_meta = row_meta))
}

# ============================== #
# 5) Run for BP / MF / CC
# ============================== #
for (ont in c("BP","MF","CC")) {
  plot_one_ontology(
    combined_long, ont_label = ont, out_prefix = out_prefix,
    min_cells = min_cells,
    top_per_parent = top_per_parent,
    top_global_cap = top_global_cap,
    value_transform_mode = value_transform_mode,
    wrap_width = wrap_width,
    row_name_max_cm = row_name_max_cm,
    cap_quantile = cap_quantile,
    cluster_cols = cluster_cols,
    show_col_dend = show_col_dend,
    family_title_fs = family_title_fs,
    pad_top_mm = pad_top_mm, pad_right_mm = pad_right_mm,
    pad_bottom_mm = pad_bottom_mm, pad_left_mm = pad_left_mm,
    left_spacer_mm = left_spacer_mm
  )
}

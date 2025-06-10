library("SingleCellExperiment")
library("scuttle")
library("edgeR")
library("purrr")
library("ggplot2")
library("EnhancedVolcano")
#library("pheatmap")
library("dplyr")
library("here")
# library("sessioninfo")

#===============================================================================
# Compute differential gene expression (DGE) test for both:
# - brain-area categorical, and
# - brain-area-numeric (linear model/GLM)  per cluster

# Overview:
# - Subsets the DGE input data (return_dge$coldata) to include only pseudo-bulk samples belonging to that cluster 
# - highlight (filter) genes that are actually expressed in each cluster accordingly with 50th quantile on edgeR::cpm(aggregated_counts)
# - On categorical var I am not computing all pairwise contrasts, only the coefficient specified
#===============================================================================
# Technical notes: scran::pseudoBulkDGE() does not allow you to specify which contrasts or groups to test directly in the function call when you have multiple coefficients (e.g., 4 brain areas), because of that, I am using manual implementation with edgeR
# scran::pseudoBulkDGE(), In short:
# It is a wrapper that helps run DGE (differential gene expression) analysis on aggregated pseudobulk data,  
# when we pass a design matrix with multiple coefficients (like 4 brain areas), it:
# - By default, return all model coefficients (not contrasts)
# - It runs a likelihood ratio test (LRT) or Wald test, depending on the method used (edgeR vs DESeq2).
#===============================================================================

#### Set up dirs ####
input_dir <- here("processed-data", "05_brain_area_differential_expression")

data_dir <- here("processed-data", "05_brain_area_differential_expression")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "05_brain_area_differential_expression", "10_pseudobulk_DEG_linear_model")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


#============= Prepare data to compute DGE with linear model trend on BayesSpace k=20 

k_merge = 20
k_merge <- sprintf("%02d", k_merge)

message("Processing BayesSpace k=", k_merge)

## Load data for the current k
data_file <- file.path(input_dir, paste0("sce_pseudo_PCA_brain_area_k", k_merge, ".rds"))
spe_data <- readRDS(data_file)

## inspect data
colnames(colData(spe_data))
# [1] "age"             "BayesSpace"      "brain_area"      "brain_area2"    
# [5] "brain_id"        "diagnosis"       "expr_chrM"       "expr_chrM_ratio"
# [9] "nspots"          "pmi"             "rin"             "sample_id"      
# [13] "sex"             "sum_umi" 
table(spe_data$BayesSpace)
# Sp20D01 Sp20D02 Sp20D03 Sp20D04 Sp20D06 Sp20D07 Sp20D08 Sp20D09 Sp20D10 Sp20D11 
#   9      12      12      12      10      12      11      12      12      11 
# Sp20D12 Sp20D13 Sp20D14 Sp20D15 Sp20D16 Sp20D17 Sp20D18 Sp20D19 Sp20D20 
#   10      11      12      12      11      12      12      11      12 
table(spe_data$brain_id, spe_data$BayesSpace)
#           Sp20D18 Sp20D19 Sp20D20
# Br8518       4       3       4
# Br9037       4       4       4
# Br9090       4       4       4
#table(spe_data$sample_id, spe_data$BayesSpace)

#===============================================================================
## merge Hb SpD(s) based on SpatialRegistration "Fine" resolution
levels(spe_data$BayesSpace)
# [1] "Sp20D01" "Sp20D02" "Sp20D03" "Sp20D04" "Sp20D06" "Sp20D07" "Sp20D08"
# [8] "Sp20D09" "Sp20D10" "Sp20D11" "Sp20D12" "Sp20D13" "Sp20D14" "Sp20D15"
# [15] "Sp20D16" "Sp20D17" "Sp20D18" "Sp20D19" "Sp20D20"

SpD20_merged <- "^Sp20D06|^Sp20D08|^Sp20D16|^Sp20D19"
spe_data$SpD20_merged <-gsub(SpD20_merged, "SpD20_Hb_merged", spe_data$BayesSpace)

# ensure all unmatched levels remain unchanged
spe_data$SpD20_merged <- ifelse(grepl(SpD20_merged, spe_data$BayesSpace),
                                "SpD20_Hb_merged", as.character(spe_data$BayesSpace))
spe_data$SpD20_merged <- factor(spe_data$SpD20_merged)
levels(spe_data$SpD20_merged)
levels(spe_data$BayesSpace)

head(table(spe_data$sample_id, spe_data$SpD20_merged))
# ...
#               SpD20_Hb_merged
# V13B23-280_A1               4
# V13B23-280_B1               4
# V13B23-280_C1               4
# V13B23-280_D1               4
# V13B23-285_A1               0
# V13B23-285_B1               3

#===============================================================================
## create new variable describing 'pseudo_brain_area' to test Habenula A-P DGE

# Br8518 = V13B23-285 defined
# Br9037 = V13B23-280
# Br9090 = V14F07-340

# add a new variable to compute DGE
assign_AP_index <- function(sample_ids) {
    case_when(
        sample_ids == "V13B23-285_A1" ~ 0,
        sample_ids %in% c("V13B23-285_B1", "V14F07-340_D1", "V13B23-280_D1") ~ 1,
        sample_ids %in% c("V13B23-285_C1", "V14F07-340_C1", "V13B23-280_C1") ~ 2,
        sample_ids %in% c("V13B23-285_D1", "V14F07-340_B1", "V14F07-340_A1", "V13B23-280_B1") ~ 3,
        sample_ids == "V13B23-280_A1" ~ 4,
        TRUE ~ NA_integer_  # fallback for unmapped samples
    )
}
# create the categorical brain-area variable
colData(spe_data)$pseudo_brain_area <- assign_AP_index(colData(spe_data)$sample_id)
# make sure all sample IDs are included
table(is.na(spe_data$pseudo_brain_area))  # should be FALSE
# duplicate variable to handle as continuous linear predictor brain-area variable
colData(spe_data)$pseudo_brain_area_numeric <- as.integer(as.character(colData(spe_data)$pseudo_brain_area))

# Make a factor for compare groups
spe_data$pseudo_brain_area <- factor(spe_data$pseudo_brain_area)
table(spe_data$BayesSpace, spe_data$pseudo_brain_area)
table(spe_data$brain_id, spe_data$pseudo_brain_area)
#         0  1  2  3  4
# Br8518 15 18 19 19  0
# Br9037  0 18 19 19 19
# Br9090  0 17 17 36  0

# creates unique factor levels for every combination of brain_id and BayesSpace and avoids accidental duplicates.
colData(spe_data)$pseudo_sample_id <- paste0(colData(spe_data)$sample_id, "_", colData(spe_data)$SpD20_merged)
# colData(spe_data)$pseudo_sample_id <- paste0(
#     colData(spe_data)$sample_id, "_", colData(spe_data)$SpD20_merged, "_", seq_len(ncol(spe_data))
# )
head(spe_data$pseudo_sample_id)
# [1] "V13B23-280_A1_Sp20D01" "V13B23-280_B1_Sp20D01" "V13B23-280_C1_Sp20D01"
# [4] "V13B23-280_D1_Sp20D01" "V13B23-285_A1_Sp20D01" "V13B23-285_B1_Sp20D01"

# inspect uniqueness
table(duplicated(colData(spe_data)$pseudo_sample_id))
# FALSE  TRUE 
# 184    32 
dups <- colData(spe_data)$pseudo_sample_id[duplicated(colData(spe_data)$pseudo_sample_id)]
# track sample_ids duplicated. Expected due H_merged is a composite region
if (length(dups)>0) {
    head(dups)
    # [1] "V13B23-280_A1_SpD20_Hb_merged" "V13B23-280_B1_SpD20_Hb_merged"
    # [3] "V13B23-280_C1_SpD20_Hb_merged" "V13B23-280_D1_SpD20_Hb_merged" ..
    dup_table <- table(colData(spe_data)$pseudo_sample_id)
    dup_table[dup_table > 1]
}

length(unique(colData(spe_data)$brain_id)) * length(unique(colData(spe_data)$SpD20_merged))
# = 3 × 16 = 48
# How many pseudo-bulk groups?
length(unique(colData(spe_data)$pseudo_sample_id))
# 216

#===============================================================================
# compute counts matrix from the global aggregation

aggregated <- scuttle::aggregateAcrossCells(
    x = spe_data,
    ids = colData(spe_data)$pseudo_sample_id,
    statistics = "sum",
    use.assay.type = "counts"
)

agg_coldata <- data.frame(
    sample_id = colnames(aggregated),  # this is the full pseudo_sample_id
    cluster = sub(".*_", "", colnames(aggregated))  # extract SpD label
)

# Now extract true sample base for mapping metadata
agg_coldata$base_id <- sub("^(([^_]+_[^_]+))_.*", "\\1", agg_coldata$sample_id)

# Add metadata (e.g. pseudo_brain_area and brain_id/donor)
map_pseudo <- as.data.frame(colData(spe_data)) |> distinct(sample_id, .keep_all = TRUE)

# add categorical brain-area variable
agg_coldata$pseudo_brain_area <- map_pseudo$pseudo_brain_area[
    match(agg_coldata$base_id, map_pseudo$sample_id)
]
# add donor to meta-data for used as covariable
agg_coldata$donor <- map_pseudo$brain_id[
    match(agg_coldata$base_id, map_pseudo$sample_id)
]
agg_coldata$donor <- factor(agg_coldata$donor)
table(agg_coldata$donor)
# Br8518 Br9037 Br9090 
# 63     63     58

#===============================================================================

plot_clusterwise_volcanos <- function(
        spe_data,
        aggregated_counts,                    # "SpatialExperiment"
        agg_coldata,                          # "data.frame"
        model_name,
        brain_area_var = "pseudo_brain_area", # or pseudo_brain_area_numeric for the linear model
        covar = "donor",                      # g.e: donor or nspots
        cluster_var = "SpD20_merged",         # relative to SpD k=20 with Hb merged clusters
        expr_threshold = 1,                   # only include genes with average CPM > 1 in the cluster (minimum expression level)
        FDR_thr = 0.05,
        output_dir = NULL,
        contrast_label,
        contrast_coef,                        # ge. coef=2 ~ area1 vs area0 
        expression_quantile = NULL
) {
    
    # expr_cpm <- edgeR::cpm(return_dge$fit$counts)
    expr_cpm <- edgeR::cpm(aggregated_counts)
    colnames(expr_cpm)
    # [1] "V13B23-280_A1_Sp20D01"         "V13B23-280_A1_Sp20D02"        
    # [3] "V13B23-280_A1_Sp20D03"         "V13B23-280_A1_Sp20D04"        
    # ... "V14F07-340_D1_SpD20_Hb_merged"
    
    # ==========================================================================
    # extract and join `nspots` from `spe_data` to `agg_coldata`
    if (covar=="nspots") {
        spot_metadata <- data.frame(
            pseudo_sample_id = colData(spe_data)$pseudo_sample_id,
            nspots = colData(spe_data)$nspots,
            stringsAsFactors = FALSE
        )
        # Ensure one row per pseudo_sample_id
        spot_metadata_unique <- spot_metadata[!duplicated(spot_metadata$pseudo_sample_id), ]
        # Join to agg_coldata by matching sample_id <-> pseudo_sample_id
        agg_coldata$nspots <- spot_metadata_unique$nspots[
            match(agg_coldata$sample_id, spot_metadata_unique$pseudo_sample_id)
        ]
        # Check for failed joins
        if (anyNA(agg_coldata$nspots)) {
            warning("Some sample_id values in agg_coldata could not be matched to nspots.")
        }
    }
    # ==========================================================================
    
    
    # Add gene names
    gene_name_map <- rowData(spe_data)$gene_name
    names(gene_name_map) <- rownames(spe_data)
    
    # Map cluster annotation
    map_cluster <- agg_coldata[, c("sample_id", "cluster")]
    head(map_cluster)
    #               sample_id cluster
    # 1 V13B23-280_A1_Sp20D01 Sp20D01
    # 2 V13B23-280_A1_Sp20D02 Sp20D02
    # 3 V13B23-280_A1_Sp20D03 Sp20D03  
    # Sanity check:
    identical(colnames(aggregated_counts), map_cluster$sample_id)
    # TRUE
    # Direct match on full pseudo sample IDs
    agg_coldata$cluster <- map_cluster$cluster[
        match(agg_coldata$sample_id, map_cluster$sample_id)
    ]
    table(agg_coldata$cluster)
    # merged Sp20D01 Sp20D02 Sp20D03 Sp20D04 Sp20D07 Sp20D09 Sp20D10 Sp20D11 Sp20D12 
    # 11       9      12      12      12      12      12      12      11      10 
    # Sp20D13 Sp20D14 Sp20D15 Sp20D17 Sp20D18 Sp20D20 
    # 11      12      12      12      12      12 
    clusters <- unique(agg_coldata$cluster)
    clusters
    # [1] "Sp20D01" "Sp20D02" "Sp20D03" "Sp20D04" "Sp20D07" "Sp20D09" "Sp20D10"
    # [8] "Sp20D11" "Sp20D12" "Sp20D13" "Sp20D14" "Sp20D15" "Sp20D17" "Sp20D18"
    # [15] "Sp20D20" "merged"
    
    message("Generating volcano plots by cluster...")
    message("Processing ", length(clusters), " levels")
    
    # pdf_file <- file.path(output_dir, paste0("volcano_", model_name, "_FDR05.pdf"))
    pdf_file <- file.path(output_dir, paste0("volcano_", model_name, "_FDR", (FDR_thr*100) ,"p_ExpQuantile", (expression_quantile*100),"th.pdf"))
    message("File name:", pdf_file)
    pdf(pdf_file, width = 8, height = 8)
    
    for (clust in clusters) {
        
        # clust="merged"
        # clust="Sp20D01"
        message("Volcano plot for cluster: ", clust)
        
        # Identify samples in cluster
        table(agg_coldata$cluster)
        samples_in_cluster <- agg_coldata$sample_id[agg_coldata$cluster == clust]
        samples_in_cluster
        # [1] "V13B23-280_A1" "V13B23-280_B1" "V13B23-280_C1" "V13B23-280_D1" ... n
        
        # Filter samples in the cluster
        sample_mask <- grepl(paste(samples_in_cluster, collapse = "|"), colnames(aggregated_counts))
        dge_cluster <- DGEList(counts = assay(aggregated_counts, "counts")[, sample_mask])
        
        # Subset agg_coldata + design, and drop levels that does not exist in that cluster to avoid error on estimateDisp()
        coldata_cluster <- droplevels(
            agg_coldata[agg_coldata$cluster == clust, ]
        )
        # # skip clusters that have fewer than 2 levels of pseudo_brain_area
        # if (nlevels(coldata_cluster$pseudo_brain_area) < 2) {
        #     message("Skipping ", clust, ": not enough brain area levels.")
        #     next
        # }
        
        # Use "brain_area_var" for categorical variable OR "brain_area_var_numeric" for linear model
        if (!is.null(covar)) {
            terms <- c(brain_area_var, covar)
        } else {
            terms <- brain_area_var
        }
        design_cluster <- model.matrix(
            reformulate(termlabels = terms),
            data = coldata_cluster
        )
        
        print(colnames(design_cluster))
        dge_cluster <- calcNormFactors(dge_cluster)
        dge_cluster <- estimateDisp(dge_cluster, design_cluster)
        fit_cluster <- glmFit(dge_cluster, design_cluster)
        lrt_cluster <- glmLRT(fit_cluster, coef = contrast_coef)  
        res_cluster <- edgeR::topTags(lrt_cluster, n = Inf)$table
        res_cluster$gene_name <- gene_name_map[rownames(res_cluster)]
        
        ## Filter for low expressed genes in the current cluster
        if (is.null(expression_quantile)) {
            cluster_samples <- coldata_cluster$sample_id
            sample_ids_cpm <- colnames(expr_cpm)
            sample_mask_cpm <- grepl(paste(cluster_samples, collapse = "|"), sample_ids_cpm)
            expressed_genes <- rownames(expr_cpm)[rowMeans(expr_cpm[, sample_mask_cpm, drop = FALSE]) > expr_threshold]
        } else {
            ## filter for low expressed genes in the current cluster using CPM quantiles
            cluster_samples <- coldata_cluster$sample_id
            sample_mask_cpm <- grepl(paste(cluster_samples, collapse = "|"), colnames(expr_cpm))
            avg_expr <- rowMeans(expr_cpm[, sample_mask_cpm, drop = FALSE])
            nonzero_expr <- avg_expr[avg_expr > 0]
            # dynamically set threshold; g.e 25th / 50th percentile threshold
            cluster_expr_threshold <- quantile(nonzero_expr, probs = expression_quantile, na.rm = TRUE)
            message("Percentile threshold set (", (expression_quantile*100), "th): ", cluster_expr_threshold)
            expressed_genes <- names(avg_expr)[avg_expr > cluster_expr_threshold]
        }        
        
        res_filtered <- res_cluster[rownames(res_cluster) %in% expressed_genes, ]
        if (nrow(res_filtered) == 0) {
            message("No expressed genes found for ", clust)
            next
        }
        
        # Identify top up/down genes based on logFC direction and significance
        #top_genes <- head(res_cluster$gene_name[order(res_cluster$FDR)], 20)
        mhb_up_genes <- res_filtered$gene_name[res_filtered$logFC > 0 & res_filtered$FDR < FDR_thr]
        lhb_up_genes <- res_filtered$gene_name[res_filtered$logFC < 0 & res_filtered$FDR < FDR_thr]
        label_genes <- c(
            head(mhb_up_genes[order(-res_filtered[res_filtered$gene_name %in% mhb_up_genes, "logFC"])], 10),
            head(lhb_up_genes[order(res_filtered[res_filtered$gene_name %in% lhb_up_genes, "logFC"])], 10)
        )
        
        # make custom subtitle
        if (!is.null(expression_quantile)) {
            sub_title <- paste0(model_name, " (FDR=", FDR_thr*100, "%; Expr.Quantile=", (expression_quantile*100) ,"th)")
        } else {
            sub_title <- paste0(model_name, " (FDR=", FDR_thr*100, "%)") 
        }
        message(sub_title)
        
        p1 <- EnhancedVolcano::EnhancedVolcano(
            res_filtered,
            lab = res_filtered$gene_name,
            selectLab = label_genes,
            x = "logFC",
            y = "FDR",
            title = paste("Habenula AP", contrast_label, "- Cluster:", clust),
            subtitle = sub_title,
            # pCutoff = 0.05, # statistical significance threshold (FDR ≤ 0.05)
            # FCcutoff = 0.5, # biological effect size threshold (log2FC > ±0.5)
            pCutoff = FDR_thr,     # even 0.2 for discovery / Default 0.05
            FCcutoff = 0.25,       # for smaller effect genes
            pointSize = 2.0,
            labSize = 4.0,
            max.overlaps = 50, 
            drawConnectors = TRUE
        )
        print(p1)
        
    }
        
    dev.off()
        
}
    
    
#===============================================================================
    
# LHb++: Column 1 (coef = 0): The intercept (baseline expression for pseudo_brain_area0)
# LHb+: Column 2 (coef = 1): log fold change of pseudo_brain_area1 relative to baseline
# LHb: Column 3 (coef = 2): logFC for pseudo_brain_area2 vs. baseline
# MHb+: Column 4 (coef = 3): logFC for pseudo_brain_area3 vs. baseline
# MHb++: Column 5 (coef = 4): logFC for pseudo_brain_area4 vs. baseline

#===============================================================================
## DGE WITH CATEGORICAL VARIABLE: Compare pseudo_brain_area0 vs pseudo_brain_area4
#===============================================================================

# From LHb+ to MHb+
for (covar_term in c("donor", "nspots")) {
    
    print(covar_term)
    model_name = paste0("BSk", k_merge, "_model1_AP_G0-G4_", covar_term)
    message("Processing categorical model: ", model_name)
    
    # build DGE and plot Volcanos
    # - FDR at 5%, with filtering by gene expr at 50% percentile, with relaxed FCcutoff for smaller effect genes
    plot_clusterwise_volcanos(
        spe_data = spe_data,
        aggregated_counts = aggregated,  # "SpatialExperiment"
        agg_coldata = agg_coldata,       # "data.frame"
        model_name = model_name,
        brain_area_var = "pseudo_brain_area",
        covar = covar_term,
        cluster_var = "SpD20_merged",
        output_dir = plot_dir,
        FDR_thr = 0.05, 
        contrast_label = paste0("G0-G4 (", covar_term, ")"),
        contrast_coef = 5, # area4 vs area0 (area0 is baseline)
        expression_quantile = 0.50 # Relaxed; allows moderately expressed genes
    )
    message("Plots done!")
    
}


# Positive logFC → genes upregulated in area4 compared to area0
# Negative logFC → genes upregulated in area0 compared to area4

#===============================================================================
## Compare pseudo_brain_area1 vs pseudo_brain_area3

# Set pseudo_brain_area1 as the reference level to compare pseudo_brain_area3 (LHb -> MHb)
# model1:  from LHb+ to MHb+
agg_coldata$pseudo_brain_area <- factor(agg_coldata$pseudo_brain_area)
agg_coldata$pseudo_brain_area <- relevel(agg_coldata$pseudo_brain_area, ref = "3")
# re-create the design matrix to check intercept
# design <- model.matrix(~ pseudo_brain_area + donor, data = agg_coldata)
# colnames(design)
# which(colnames(design) == "pseudo_brain_area1")

for (covar_term in c("donor", "nspots")) {
    
    print(covar_term)
    model_name = paste0("BSk", k_merge,  "_model2_AP_G1-G3_", covar_term)
    message("Processing categorical model: ", model_name)
    
    # build DGE and plot Volcanos
    # - FDR at 5%, with filtering by gene expr at 50% percentile, with relaxed FCcutoff for smaller effect genes
    plot_clusterwise_volcanos(
        spe_data = spe_data,
        aggregated_counts = aggregated, #  assay(aggregated, "counts"),
        agg_coldata = agg_coldata,
        model_name = model_name,
        brain_area_var = "pseudo_brain_area",
        covar = covar_term,
        cluster_var = "SpD20_merged",
        output_dir = plot_dir,
        FDR_thr = 0.05,            # try 0.1 for EDA or event 0.2 for discovery
        contrast_label = paste0("G1-G3 (", covar_term, ")"),
        contrast_coef = 2, 
        expression_quantile = 0.50 # 0.50 Relaxed FOR EDA;  0.75 a bit more strict; 0.97 most to ideal=1 
    )
    message("Plots done!")
    
}


#===============================================================================
# DGE LINEAR MODEL: pseudo_brain_area is modeled as a linear numeric variable,
# enabling tests of trend along the AP axis (e.g., increasing index)
#===============================================================================

## add to pseudo_brain_area_numeric to meta-data
agg_coldata$pseudo_brain_area_numeric <- assign_AP_index(agg_coldata$base_id)
# make sure all sample IDs are included
table(is.na(spe_data$pseudo_brain_area_numeric))  # should be FALSE

for (covar_term in c("donor", "nspots")) {
    
    print(covar_term)
    model_name = paste0("BSk", k_merge,  "_model3_linear_AP_0-4_", covar_term)
    message("Processing linear model: ", model_name)
    
    # build DGE and plot Volcanos
    # - FDR at 5%, with filtering by gene expr at 50th percentile, with relaxed FCcutoff for smaller effect genes
    plot_clusterwise_volcanos(
        spe_data = spe_data,
        aggregated_counts = aggregated, 
        agg_coldata = agg_coldata,
        model_name = model_name,
        brain_area_var = "pseudo_brain_area_numeric", 
        covar = covar_term,
        cluster_var = "SpD20_merged",
        output_dir = plot_dir,
        FDR_thr = 0.05,
        contrast_label = paste0("linear:0-4 (", covar_term, ")"),
        contrast_coef = 2,
        expression_quantile = 0.50
    )
    message("Plots done!")
    
}


#===============================================================================
# compute DGE and plot volcano plots only for samples with pseudo_brain_area_numeric from 1 to 3

# Subset only pseudo_brain_area_numeric in 1, 2, 3
subset_ids <- agg_coldata$sample_id[agg_coldata$pseudo_brain_area_numeric %in% 1:3]
agg_coldata_subset <- agg_coldata[agg_coldata$sample_id %in% subset_ids, ]
# Subset counts
aggregated_subset <- aggregated[, colnames(aggregated) %in% subset_ids]
# Subset the main SPE object if needed (for rowData)
spe_data_subset <- spe_data[, colnames(spe_data) %in% subset_ids]

for (covar_term in c("donor", "nspots")) {
    
    print(covar_term)
    model_name = paste0("BSk", k_merge, "_model4_linear_AP_1-3_", covar_term)
    message("Processing linear model: ", model_name)
    
    # build DGE and plot Volcanos
    # - FDR at 5%, with filtering by gene expr at 50th percentile, with relaxed FCcutoff for smaller effect genes
    plot_clusterwise_volcanos(
        spe_data = spe_data_subset,
        aggregated_counts = aggregated_subset,
        agg_coldata = agg_coldata_subset,
        model_name = model_name,
        brain_area_var = "pseudo_brain_area_numeric", 
        covar = covar_term,
        cluster_var = "SpD20_merged",
        output_dir = plot_dir,
        FDR_thr = 0.05,
        contrast_label = paste0("linear:1-3 (", covar_term, ")"),
        contrast_coef = 2,
        expression_quantile = 0.50
    )
    message("Plots done!")
    
}

# Need to check: nspots in linear AP 1-3 Fails
# Error in glmFit.default(sely, design, offset = seloffset, dispersion = 0.05,  : 
#                             nrow(design) disagrees with ncol(y)
#                         In addition: Warning message:
#                             In plot_clusterwise_volcanos(spe_data = spe_data_subset, aggregated_counts = aggregated_subset,  :
#                                                              Some sample_id values in agg_coldata could not be matched to nspots.



# model_name = paste0("BSk", k_merge, "_model4_linear_AP_1-3_nspots")
# model_name
# ## check, fails
# plot_clusterwise_volcanos(
#     spe_data = spe_data_subset,
#     aggregated_counts = aggregated_subset,
#     agg_coldata = agg_coldata_subset,
#     model_name = model_name,
#     brain_area_var = "pseudo_brain_area_numeric", 
#     covar = "nspots", 
#     cluster_var = "SpD20_merged",
#     output_dir = plot_dir,
#     FDR_thr = 0.05,
#     contrast_label = "linear (donor)",
#     contrast_coef = 2,
#     expression_quantile = 0.50
# )


    
    #===============================================================================
    
    # top_genes <- rownames(topTags(return_dge$lrt, n = 50)$table)
    # dge <- return_dge$fit
    # counts <- cpm(dge$counts, log = TRUE)  # log2 CPM
    # 
    # # Subset to top genes
    # heatmap_matrix <- counts[top_genes, ]
    # # Z-score normalize by gene (row-wise)
    # heatmap_matrix_z <- t(scale(t(heatmap_matrix)))
    # 
    # # Add annotation for brain_area2_numeric (for columns)
    # sample_metadata <- return_dge$coldata
    # rownames(sample_metadata) <- colnames(heatmap_matrix)
    # 
    # annotation_col <- data.frame(
    #     brain_area2_numeric = sample_metadata$pseudo_brain_area
    # )
    # rownames(annotation_col) <- colnames(heatmap_matrix)
    # 
    # 
    # ## prepare heatmaps with gene-expr data
    # 
    # pheatmap(
    #     heatmap_matrix_z,
    #     annotation_col = annotation_col,
    #     cluster_rows = TRUE,
    #     cluster_cols = TRUE,
    #     show_rownames = TRUE,
    #     show_colnames = FALSE,
    #     fontsize_row = 6,
    #     main = "Top DE Genes Heatmap"
    # )
    # 
    # 
    # # Top 50 genes (already from previous steps)
    # top_genes <- rownames(edgeR::topTags(return_dge$lrt, n = 50)$table)
    # 
    # # LogCPM from DGEList
    # dge <- return_dge$fit
    # logCPM <- edgeR::cpm(dge$counts, log = TRUE)
    # 
    # # Subset to top genes
    # heatmap_matrix <- logCPM[top_genes, ]
    # 
    # # Z-score normalization (gene-wise)
    # heatmap_matrix_z <- t(scale(t(heatmap_matrix)))
    # 
    # # Column annotation
    # sample_metadata <- de_results$coldata
    # rownames(sample_metadata) <- colnames(heatmap_matrix_z)
    # 
    # annotation_col <- data.frame(
    #     pseudo_brain_area = sample_metadata$pseudo_brain_area
    # )
    # rownames(annotation_col) <- colnames(heatmap_matrix_z)
    # 
    # # Plot heatmap with clustering enabled
    # pheatmap(
    #     heatmap_matrix_z,
    #     annotation_col = annotation_col,
    #     cluster_rows = TRUE,     # cluster genes
    #     cluster_cols = TRUE,     # cluster samples
    #     show_rownames = TRUE,
    #     show_colnames = FALSE,
    #     cutree_cols = 3,
    #     fontsize_row = 6,
    #     main = "Top 50 DE Genes (Clustered)"
    # )
    
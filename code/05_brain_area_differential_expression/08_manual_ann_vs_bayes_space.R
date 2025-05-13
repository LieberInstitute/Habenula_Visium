library("spatialLIBD")
library("dplyr")
library("ggplot2")
library("gridExtra")
library("readr")
library("here")
library("sessioninfo")


k <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# k = 13
k_nice <- sprintf("%02d", k) # Format k
message("Processing BayesSpace k=", k_nice)

dir_plots <- here("plots", "05_brain_area_differential_expression", "08_manual_ann_vs_bayes_space")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

manual_ann_dir <- here("processed-data", "03_spatialLIBD_app", "Manual_annotations")
rds_dir <- here("processed-data", "04_harmony_BayesSpace", "spe_harmony.rds")
spe <- readRDS(rds_dir)

## Import BayesSpace clusters
clusters_BayesSpace_dir <- here(
    "processed-data",
    "04_harmony_BayesSpace",
    "clusters_BayesSpace"
)
# head(spe$key[1:5])
spe <- cluster_import(
    spe,
    cluster_dir = clusters_BayesSpace_dir,
    prefix = "",
    overwrite = TRUE
)
# Overwriting 'spe$key'. Set 'overwrite = FALSE' if you do not want to overwrite it.


## quick inspection: check how many genes expressed by cluster we have after pseudobulk

spe
colnames(colData(spe))
spot_names <- colnames(spe)
head(spot_names)
# [1] "AAACAAGTATCTCCCA-1" "AAACACCAATAACTGC-1" "AAACAGCTTTCAGAAG-1"
# [4] "AAACAGGGTCTATATT-1" "AAACAGTGTTCCTGGG-1" "AAACATGGTGAGAGGA-1"

## prepare manual annotation

f_name <- here(manual_ann_dir, "spatialLIBD_ManualAnnotation_2025-02-07_KDM-CSC_RNAScope_merged_spatialLIBD.csv")
#readLines(f_name)
df_manual_annotations <- as.data.frame(read_delim(f_name, delim = ",", show_col_types = FALSE))
dim(df_manual_annotations)
# [1] 2089    3
head(df_manual_annotations)                                                                                      
# sample_id          spot_name          ManualAnnotation
# 1 V13B23-285_A1 AAACGAAGAACATACC-1      Broad_Hb_v1
# 2 V13B23-285_A1 AAATAACCATACGGGA-1      Broad_Hb_v1
# 3 V13B23-285_A1 AAATGGCATGTCTTGT-1      Broad_Hb_v1
# 4 V13B23-285_A1 AAATGTGGGTGCTCCT-1      Broad_Hb_v1
# 5 V13B23-285_A1 AAATTACACGACTCTG-1      Broad_Hb_v1
# 6 V13B23-285_A1 AACACGACTGTACTGA-1      Broad_Hb_v1    
unique(df_manual_annotations$ManualAnnotation)

message("Computing Hb domains for BS k: ", hab_level)
    
## extract meta-data
# colnames(spe)
# colData(spe)
df_domains <- colData(spe)
head(df_domains)
dim(df_domains)

df_domains$spot_name_ann <- row.names(df_domains)
row.names(df_domains) <- NULL
df_domains <- as.data.frame(df_domains)
colnames(df_domains)

## extract BayesSpace_harmony names
all_domains <- colnames(df_domains)[grep("BayesSpace_harmony", colnames(df_domains))]
all_domains

# Set PDF for combine plot by sample and cluster in x-axis
pdf(file = file.path(dir_plots, paste0("stacked_bar_manual-vs-BS_by_clusters.pdf")))

for (SpD in all_domains) {
    
    # SpD = "BayesSpace_harmony_k03"
    message("Searching matching spots between manual annotations and ", SpD)
    
    # select columns to use
    df_domain <- df_domains |>
        select(all_of(c("sample_id", SpD, "brain_id", "spot_name_ann")))
    print(head(df_domain, n=3))
    #     sample_id BayesSpace_harmony_k13 brain_id      spot_name_ann
    # 1 V13B23-280_A1                      8   Br9037 AAACAAGTATCTCCCA-1
    # 2 V13B23-280_A1                     11   Br9037 AAACACCAATAACTGC-1
    # 3 V13B23-280_A1                      8   Br9037 AAACAGCTTTCAGAAG-1
    message("Total spots: ", length(df_domain$spot_name_ann))
    
    # Mark Matches and Non-Matches by sample
    df_domain_labeled <- df_domain %>%
        mutate(match_status = ifelse(spot_name_ann %in% df_manual_annotations$spot_name, "Habenula", "No-Habenula")) %>%
        left_join(
            df_manual_annotations[, c("spot_name", "ManualAnnotation")],
            by = c("spot_name_ann" = "spot_name"),
            relationship = "many-to-many"
        )
    #names(df_domain_labeled)
    #head(df_domain_labeled)
    #table(df_domain_labeled$match_status, df_domain_labeled$sample_id)
    
    # count of how many spots per ManualAnnotation are Habenula or No-Habenula
    df_plot <- df_domain_labeled %>%
        group_by(brain_id, sample_id, !!sym(SpD), match_status) %>%
        summarise(count = n(), .groups = "drop") # plot absolute counts
    # summarise(prop = n() / sum(n()), .groups = "drop") # plot proportions
    head(df_plot)

    # # Plot #spots by brain_id
    # ggplot(df_plot, aes(x = as.factor(!!sym(SpD)), y = count, fill = match_status)) +
    #     geom_bar(stat = "identity", position = "stack") +
    #     facet_wrap(~ brain_id) +
    #     labs(
    #         title = paste0("Habenula vs No-Habenula by ", SpD, " and Sample"),
    #         x = SpD,
    #         y = "Number of Spots",
    #         fill = "Match Status"
    #     ) +
    #     theme_minimal() +
    #     theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Combine sample and cluster in x-axis
    plt1 <- ggplot(df_plot, aes(x = interaction(sample_id, !!sym(SpD)), y = count, fill = match_status)) +
        geom_bar(stat = "identity") +
        labs(
            title = paste0("Habenula vs No-Habenula by ", SpD, " and Sample"),
            x = "Sample + Cluster",
            y = "Number of Spots",
            fill = "Match Status"
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    print(plt1)
    
}

dev.off()    


# Set PDF for combine plot by sample and cluster in x-axis
pdf(file = file.path(dir_plots, paste0("stacked_bar_manual-vs-BS_by_clusters.pdf")))

for (SpD in all_domains) {
    
    # SpD = "BayesSpace_harmony_k03"
    message("Searching matching spots between manual annotations and ", SpD)
    
    # select columns to use
    df_domain <- df_domains |>
        select(all_of(c("sample_id", SpD, "brain_id", "spot_name_ann")))
    print(head(df_domain, n=3))
    #     sample_id BayesSpace_harmony_k13 brain_id      spot_name_ann
    # 1 V13B23-280_A1                      8   Br9037 AAACAAGTATCTCCCA-1
    # 2 V13B23-280_A1                     11   Br9037 AAACACCAATAACTGC-1
    # 3 V13B23-280_A1                      8   Br9037 AAACAGCTTTCAGAAG-1
    message("Total spots: ", length(df_domain$spot_name_ann))
    
    # Mark Matches and Non-Matches by sample
    df_domain_labeled <- df_domain %>%
        mutate(match_status = ifelse(spot_name_ann %in% df_manual_annotations$spot_name, "Habenula", "No-Habenula")) %>%
        left_join(
            df_manual_annotations[, c("spot_name", "ManualAnnotation")],
            by = c("spot_name_ann" = "spot_name"),
            relationship = "many-to-many"
        )

    # compute proportions
    df_plot <- df_domain_labeled %>%
        group_by(brain_id, sample_id, !!sym(SpD), match_status) %>%
        summarise(count = n(), .groups = "drop") %>%
        group_by(brain_id, sample_id, !!sym(SpD)) %>%
        mutate(prop = count / sum(count)) %>%
        ungroup()
    #head(df_plot)
    
    # plot stacked bar with proportions by cluster in the x-axis
    plt2 <- ggplot(df_plot, aes(x = as.factor(.data[[SpD]]), y = prop, fill = match_status)) +
        geom_bar(stat = "identity") +
        labs(
            title = paste0("Proportion of Habenula vs No-Habenula by ", SpD),
            x = SpD,
            y = "Proportion of Spots",
            fill = "Match Status"
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+
    # geom_text(
    #     aes(label = scales::percent(prop, accuracy = 1)),
    #     position = position_stack(vjust = 0.5),
    #     size = 3,
    #     color = "black"
    # )
    
    
    print(plt2)
    
}

dev.off()    


message(' Plots completed!')





## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
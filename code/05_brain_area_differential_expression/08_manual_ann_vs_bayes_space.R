########################################################################
## Build stacked bar plots to compare Habenula manual annotations against BayesSpace clusters
## Authors. CSC
## Data: May 13rd, 2025
## For 60 to 80k spots: $srun --pty --mem=30GB --x11 bash
########################################################################


library("spatialLIBD")
library("dplyr")
library("stringr")
library("ggplot2")
library("gridExtra")
library("readr")
library("here")
library("sessioninfo")


dir_plots <- here("plots", "05_brain_area_differential_expression", "08_manual_ann_vs_bayes_space")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

manual_ann_dir <- here("processed-data", "03_spatialLIBD_app", "Manual_annotations")
rds_dir <- here("processed-data", "05_brain_area_differential_expression", "spe_harmony_ann.rds")

spe <- readRDS(rds_dir)
spe

## quick inspection: check how many genes expressed by cluster we have after pseudobulk

colnames(colData(spe))
levels(colData(spe)$BayesSpace)
colnames(colData(spe))[grep("BayesSpace_harmony_", colnames(colData(spe)))]

# first make manual verification to identify missing SpD(s) in samples causing visualization issues
table(colData(spe)$sample_id, colData(spe)$BayesSpace)
spot_names <- colnames(spe)
head(spot_names)
# [1] "AAACAAGTATCTCCCA-1" "AAACACCAATAACTGC-1" "AAACAGCTTTCAGAAG-1"
# [4] "AAACAGGGTCTATATT-1" "AAACAGTGTTCCTGGG-1" "AAACATGGTGAGAGGA-1"

## prepare Hb RNAScope Manual Annotations (KDM)

f_name <- here(manual_ann_dir, "spatialLIBD_ManualAnnotation_2025-02-07_KDM-CSC_RNAScope_merged_spatialLIBD.csv")
#readLines(f_name)
df_manual_annotations <- as.data.frame(read_delim(f_name, delim = ",", show_col_types = FALSE))

## delete the spots annotated no longer use it. eg. bad quality samples
length(df_manual_annotations$sample_id)
unique(df_manual_annotations$sample_id)
df_manual_annotations <- df_manual_annotations[!grepl("V13B23-281_[A-D]1", df_manual_annotations$sample_id), ]
unique(df_manual_annotations$sample_id)

# add a unique column to match barcodes
df_manual_annotations$spot_name2 <- paste0(df_manual_annotations$sample_id, "-", df_manual_annotations$spot_name)
length(df_manual_annotations$ManualAnnotation)
message("Spots manually annotated: ", length(df_manual_annotations$ManualAnnotation))
# Spots manually annotated: 1190

#head(df_manual_annotations)                                                                                      

## extract meta-data
df_domains <- colData(spe)
#head(df_domains)
dim(df_domains)

df_domains$spot_name_ann <- row.names(df_domains)
row.names(df_domains) <- NULL
df_domains <- as.data.frame(df_domains)
#colnames(df_domains)

## extract BayesSpace_harmony names
all_domains <- colnames(df_domains)[grep("BayesSpace_harmony", colnames(df_domains))]
#all_domains

# Set PDF for combine plot by sample and cluster in x-axis
pdf(file = file.path(dir_plots, paste0("histogram_bar_manual-vs-BS_by_clusters-sample.pdf")))

for (SpD in all_domains) {
    
    # SpD = all_domains[5]
    message("Searching matching spots between manual annotations and ", SpD)
    # extract BS-k to separate in the plots by bins
    k <- str_extract(SpD, "\\d+")
    k <- as.numeric(k)
    
    # select columns to use
    df_domain <- df_domains |>
        select(all_of(c("sample_id", SpD, "brain_id", "spot_name_ann")))
    
    # add a unique column to match barcodes
    df_domain$spot_name_ann2 <- paste0(df_domain$sample_id, "-", df_domain$spot_name_ann)
    #print(head(df_domain, n=3))
    # sample_id BayesSpace_harmony_k06 brain_id      spot_name_ann
    # 1 V13B23-280_A1                      3   Br9037 AAACAAGTATCTCCCA-1
    # 2 V13B23-280_A1                      3   Br9037 AAACACCAATAACTGC-1
    # 3 V13B23-280_A1                      3   Br9037 AAACAGCTTTCAGAAG-1
    # spot_name_ann2
    # 1 V13B23-280_A1-AAACAAGTATCTCCCA-1
    # 2 V13B23-280_A1-AAACACCAATAACTGC-1
    # 3 V13B23-280_A1-AAACAGCTTTCAGAAG-1
    
    # remove samples not manually annotated - to match with the samples annotated
    if (!length(unique(df_domain$sample_id)) == length(unique(df_manual_annotations$sample_id))) {
        df_domain <- df_domain |>
            filter(sample_id %in% unique(df_manual_annotations$sample_id))
    }    
    unique(df_domain$sample_id)
    unique(df_manual_annotations$sample_id)
    message("Total spots: ", length(df_domain$spot_name_ann))
    
    # Mark Matches and Non-Matches by sample
    anyDuplicated(df_domain$spot_name_ann2)         # should be 0
    anyDuplicated(df_manual_annotations$spot_name2) # should be 0
    #head(df_domain)
    df_domain_labeled <- df_domain |>
        mutate(match_status = ifelse(spot_name_ann2 %in% df_manual_annotations$spot_name2, "Habenula", "No-Habenula")) |>
        left_join(
            df_manual_annotations[, "spot_name2", drop = FALSE],
            by = c("spot_name_ann2" = "spot_name2")
            #relationship = "many-to-many" -> spots are only unique by capture area
        )
    #names(df_domain_labeled)
    #head(df_domain_labeled)
    #table(df_domain_labeled$match_status, df_domain_labeled$sample_id)
    
    # count of how many spots per ManualAnnotation are Habenula or No-Habenula
    df_plot <- df_domain_labeled |>
        group_by(brain_id, sample_id, !!sym(SpD), match_status) |>
        summarise(count = n(), .groups = "drop") # plot absolute counts
    #head(df_plot)
    
    # Combine sample and cluster in x-axis
    df_plot <- df_plot |>
        mutate(x_label = interaction(sample_id, !!sym(SpD)))

    # Extract cluster values (e.g., 1–13) from the SpD column
    cluster_vals <- df_plot[[SpD]]
    x_labels <- levels(df_plot$x_label)
    
    # Build label vector: show cluster value every k bars
    custom_labels <- ifelse(seq_along(x_labels) %% k == 1, as.character(cluster_vals), "")
    
    # Get cluster values in the same order as x_label levels
    cluster_by_x <- df_plot |>
        distinct(x_label, .keep_all = TRUE) |>
        arrange(factor(x_label, levels = levels(df_plot$x_label))) |>
        pull(!!sym(SpD))
    
    # prepare customized x-axis labels and add vertical lines to separate clusters
    # Keep the first occurrence of each new cluster, blank for repeated values
    custom_labels <- ifelse(
        c(TRUE, diff(cluster_by_x) != 0),  # first is always TRUE
        as.character(cluster_by_x),
        ""
    )
    custom_labels <- ifelse(custom_labels != "", paste0("SpD ", custom_labels), "")
    
    # Define bar positions based dynamically
    bar_positions <- which(custom_labels != "")-1

    plt1 <- ggplot(df_plot, aes(x = x_label, y = count, fill = match_status)) +
        geom_bar(stat = "identity") +
        geom_vline(xintercept = bar_positions + 0.5, linetype = "solid", color = "darkgray") +
        scale_x_discrete(labels = custom_labels) +
        labs(
            title = paste0("Habenula vs No-Habenula: ", SpD),
            x = paste0("Spatial-Domains"),
            y = "Number of Spots",
            fill = "Match Status"
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            legend.position = "top"
        )
    print(plt1)
    
        
}

dev.off() 

message("Histogram done!")


# Set PDF for combine plot by sample and cluster in x-axis
pdf(file = file.path(dir_plots, paste0("stacked_bar_manual-vs-BS_by_clusters.pdf")))

for (SpD in all_domains) {
    
    # SpD = all_domains[10]
    message("Searching matching spots between manual annotations and ", SpD)
    
    # select columns to use
    df_domain <- df_domains |>
        select(all_of(c("sample_id", SpD, "brain_id", "spot_name_ann")))
    # add a unique column to match barcodes
    df_domain$spot_name_ann2 <- paste0(df_domain$sample_id, "-", df_domain$spot_name_ann)
    #print(head(df_domain, n=3))

    message("Total spots: ", length(df_domain$spot_name_ann))
    
    # Mark Matches and Non-Matches by sample
    df_domain_labeled <- df_domain|>
        mutate(match_status = ifelse(spot_name_ann2 %in% df_manual_annotations$spot_name2, "Habenula", "No-Habenula"))
    
    df_plot <- df_domain_labeled |>
        group_by(brain_id, sample_id, !!sym(SpD), match_status) |>
        summarise(count = n(), .groups = "drop") |>
        group_by(brain_id, sample_id, !!sym(SpD)) |>
        mutate(prop = count / sum(count)) |>
        ungroup()

    # aggregate
    df_plot2 <- df_plot |>
        group_by(!!sym(SpD), match_status) |>
        summarise(prop = mean(prop), .groups = "drop") |>
        ungroup()
    
    df_plot2 <- df_plot2 |>
        mutate(
            SpD_label = paste0("SpD", stringr::str_pad(.data[[SpD]], width = 2, pad = "0"))
        )
    
    # plot stacked bar with proportions by cluster in the x-axis
    plt1 <- ggplot(df_plot2, aes(x = SpD_label, y = prop, fill = match_status)) +
        geom_bar(position = "fill", stat = "identity") +
        scale_y_continuous(
            breaks = seq(0, 1, by = 0.25),   
            labels = scales::number_format(accuracy = 0.01)
        ) +
        labs(
            title = paste0("Habenula vs No-Habenula: ", SpD),
            x = NULL,
            y = "Proportion of Hb Spots",
            fill = "Match Status"
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(size = 12, angle = 45),
            axis.text.y = element_text(size = 12),
            plot.title   = element_text(size = 14, face = "bold", colour = "black"),
            legend.title = element_text(colour = "black"),
            legend.text  = element_text(size = 12, colour = "black"),
            legend.position = "bottom"
        ) 
    print(plt1)
    
}

dev.off()   

message("StackedPlot done!")


## =============================================================================
## Additional analysis for BS k=28

k_list <- c(11, 15, 20, 24, 28)

for (k in k_list) {

    ## select BayesSpace_harmony of interest
    BS_k <- paste0("BayesSpace_harmony_k", k)
    
    message("Processing histogram for Hb SpD(s) ", BS_k)
    
    df_domains_k <- df_domains[, c("sample_id", "brain_id", "spot_name_ann", BS_k)]
    #df_domains_k
    # Keep rows where BayesSpace_harmony_k is in the vector c(5, 10, 20, 27)
    unique(df_domains_k[[BS_k]])
    
    # Define your filtering condition based on 'k'
    df_domain_filtered <- switch(as.character(k),
                                 "11" = subset(df_domains_k, BayesSpace_harmony_k11 %in% c(6, 11)),
                                 "15" = subset(df_domains_k, BayesSpace_harmony_k15 %in% c(6, 10, 14)),
                                 "20" = subset(df_domains_k, BayesSpace_harmony_k20 %in% c(6, 8, 16, 19)),
                                 "24" = subset(df_domains_k, BayesSpace_harmony_k24 %in% c(7, 18, 22, 23, 24)),
                                 "28" = subset(df_domains_k, BayesSpace_harmony_k28 %in% c(5, 10, 11, 20, 27)),
                                 stop("SpD(s) not specified")
    )
        
    # add a unique column to match barcodes
    df_domain_filtered$spot_name_ann2 <- paste0(df_domain_filtered$sample_id, "-", df_domain_filtered$spot_name_ann)
    head(df_domain_filtered, n=3)
        
    # remove samples not manually annotated - to match with the samples annotated
    if (!length(unique(df_domain_filtered$sample_id)) == length(unique(df_manual_annotations$sample_id))) {
        df_domain_filtered <- df_domain_filtered |>
            filter(sample_id %in% unique(df_manual_annotations$sample_id))
    }    
    unique(df_domain_filtered$sample_id)
    unique(df_manual_annotations$sample_id)
    
    spots_total <- length(df_domain_filtered$spot_name_ann)
    message("Total spots: ", spots_total)
    
    # Mark Matches and Non-Matches by sample
    anyDuplicated(df_domain_filtered$spot_name_ann2)         # should be 0
    anyDuplicated(df_manual_annotations$spot_name2) # should be 0
    
    df_domain_labeled <- df_domain_filtered |>
        mutate(match_status = ifelse(spot_name_ann2 %in% df_manual_annotations$spot_name2, "Habenula", "No-Habenula")) |>
        left_join(
            df_manual_annotations[, "spot_name2", drop = FALSE],
            by = c("spot_name_ann2" = "spot_name2")
        )
    #names(df_domain_labeled)
    #head(df_domain_labeled)
    #table(df_domain_labeled$match_status, df_domain_labeled$sample_id)
    
    # count of how many spots per ManualAnnotation are Habenula or No-Habenula
    SpD <- BS_k
    df_plot <- df_domain_labeled |>
        group_by(brain_id, sample_id, !!sym(SpD), match_status) |>
        summarise(count = n(), .groups = "drop") # plot absolute counts
    #head(df_plot)
    ## calculate number of Habenula spots vs not-habneula (to use on subtitle plots)
    df_summary <- df_plot |>
        group_by(match_status) |>
        summarise(total_count = sum(count), .groups = "drop")
    
    # Combine sample and cluster in x-axis
    df_plot <- df_plot |>
        mutate(x_label = interaction(sample_id, !!sym(SpD)))
    
    # Extract cluster values (e.g., 1–13) from the SpD column
    cluster_vals <- df_plot[[SpD]]
    x_labels <- levels(df_plot$x_label)
    
    # Build label vector: show cluster value every k bars
    custom_labels <- ifelse(seq_along(x_labels) %% k == 1, as.character(cluster_vals), "")
    
    # Get cluster values in the same order as x_label levels
    cluster_by_x <- df_plot |>
        distinct(x_label, .keep_all = TRUE) |>
        arrange(factor(x_label, levels = levels(df_plot$x_label))) |>
        pull(!!sym(SpD))
    
    # prepare customized x-axis labels and add vertical lines to separate clusters
    # Keep the first occurrence of each new cluster, blank for repeated values
    custom_labels <- ifelse(
        c(TRUE, diff(cluster_by_x) != 0),  # first is always TRUE
        as.character(cluster_by_x),
        ""
    )
    custom_labels <- ifelse(custom_labels != "", paste0("SpD ", custom_labels), "")
    
    # Include sample_id in the x_label labels for the plot
    
    custom_labels <- paste0(custom_labels, " (", gsub(":.*", "", sub("\\..*", "", x_labels)), ")")
    
    # Define bar positions based dynamically
    bar_positions <- which(custom_labels != "")-1
    
    # Define custom subtitle
    spots_no_hb <- paste("No-Habenula counts =", df_summary$total_count[df_summary$match_status=="No-Habenula"])
    spots_hb <- paste(spots_no_hb, "\nHabenula counts =", df_summary$total_count[df_summary$match_status=="Habenula"])
    
    plt1 <- ggplot(df_plot, aes(x = x_label, y = count, fill = match_status)) +
        geom_bar(stat = "identity") +
        geom_vline(xintercept = bar_positions + 0.5, linetype = "solid", color = "darkgray") +
        scale_x_discrete(labels = custom_labels) +
        labs(
            title = paste0("Habenula vs No-Habenula: ", SpD),
            subtitle = paste0(spots_hb), 
            x = paste0("Spatial-Domains"),
            y = "Number of Spots",
            fill = "Match Status"
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, size=8, hjust = 1),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            legend.position = "top"
        )
    
    # Set PDF for combine plot by sample and cluster in x-axis
    pdf(file = file.path(dir_plots, paste0("Habenula_BS_k", k, "_histogram_bar_manual-vs-BS_by_clusters-sample.pdf")))
    #print(plt1)
    dev.off() 
    
    message("Habenula Histogram done!")

}

## =============================================================================




## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
library("spatialLIBD")
library("SpatialExperiment")
library("here")
library("tidyverse")
library("sessioninfo")



#######For EDA plot gene-distribution of Visium Samples #######

dir_RDS <- here("processed-data", "02_build_spe", "spe_scran_spotsweeper.rds")
dir_plots <- here("plots", "02_build_spe")

set.seed(20250221)

## load a not-filtered spe object
spe <- readRDS(dir_RDS)

## Verified number of TRUE spots in tissue
#in_tissue_spots <- sum(as.numeric(map(unique(spe$sample_id), ~ sum(spe$in_tissue[spe$sample_id == .x]))))
message("Spots in tissue: ", length(spe$in_tissue))
# Spots in tissue: 33587

# Sum UMI counts per spot
sum_umi <- colSums(counts(spe))
# Extract sample information from 'sample_id' column in colData
sample_info <- colData(spe)$sample_id
# Create a data frame
df <- data.frame(SumUMI = sum_umi, Sample = as.factor(sample_info))
head(df)

## Plot gene distribution histogram and density based on UMIs by sample

pdf(file.path(dir_plots, "distribution_sumUMI_sample.pdf"), useDingbats = FALSE)

# As distribution is highly skewed, we use a log transformation histogram
ggplot(df, aes(x = log10(SumUMI + 1), fill = Sample)) +
  geom_histogram(binwidth = 0.1, alpha = 0.6, position = "identity") +
  labs(title = "Log-transformed Distribution of Sum UMI Counts per Spot by Sample",
       x = "Log10(Sum UMI Counts + 1)",
       y = "Frequency") +
  theme_bw() +
  theme(legend.position="none") +
  # theme(axis.text = element_text(size = 6)) +
  facet_wrap(~ Sample, scales = "free_y")


## Log-transformed UMI distribution in one density plot

# Create a data frame
df <- data.frame(LogSumUMI = log10(sum_umi + 1), Sample = as.factor(sample_info))

## Plot density distribution with distinct colors per sample
ggplot(df, aes(x = LogSumUMI, color = Sample, fill = Sample)) +
  geom_density(alpha = 0.3) +  # Semi-transparent fill for visibility
  scale_fill_brewer(palette = "Set2") +  # Nice color palette for distinction
  scale_color_brewer(palette = "Set2") +
  labs(title = "Log-transformed UMI Count Distribution by Sample",
       x = "Log10(Sum UMI Counts + 1)",
       y = "Density") +
  theme_bw() +
  theme(legend.position = "bottom")  



## Box-plot with jitter points by sample 
ggplot(df, aes(x = Sample, y = LogSumUMI, fill = Sample)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +  # Boxplot without outliers for clarity
  geom_jitter(shape = 21, size = 0.5, alpha = 0.5, width = 0.2) +  # Add jittered points
  scale_fill_brewer(palette = "Set2") +  # Nice color palette
  labs(title = "Log-Transformed Sum UMI Counts per Sample",
       x = "Sample",
       y = "Log10(Sum UMI Counts + 1)") +
  theme_bw() +
  theme(axis.text = element_text(size = 8, angle = 90)) +
  theme(legend.position = "none") 


dev.off()


## Plot gene distribution histogram and density based on gene-sum by sample

# Calculate the number of detected genes per spot (genes with nonzero counts)

colnames(colData(spe))
gene_counts_per_spot <- spe$sum_gene
# Create a data frame
df <- data.frame(SumGene = gene_counts_per_spot, Sample = as.factor(sample_info))
head(df)

pdf(file.path(dir_plots, "distribution_sumGene_sample.pdf"), useDingbats = FALSE)

# As distribution is highly skewed, we use a log transformation histogram
ggplot(df, aes(x = log10(SumGene + 1), fill = Sample)) +
  geom_histogram(binwidth = 0.1, alpha = 0.6, position = "identity") +
  labs(title = "Log-transformed Distribution of Sum Gene Counts per Spot by Sample",
       x = "Log10(Sum Gene Counts + 1)",
       y = "Frequency") +
  theme_bw() +
  theme(legend.position="none") +
  # theme(axis.text = element_text(size = 6)) +
  facet_wrap(~ Sample, scales = "free_y")


## Log-transformed gene distribution in one density plot

# Create a data frame
df <- data.frame(LogSumGene = log10(gene_counts_per_spot + 1), Sample = as.factor(sample_info))

## Plot density distribution with distinct colors per sample
ggplot(df, aes(x = LogSumGene, color = Sample, fill = Sample)) +
  geom_density(alpha = 0.3) +  # Semi-transparent fill for visibility
  scale_fill_brewer(palette = "Set2") +  # Nice color palette for distinction
  scale_color_brewer(palette = "Set2") +
  labs(title = "Log-transformed Gene Count Distribution by Sample",
       x = "Log10(Sum Gene Counts + 1)",
       y = "Density") +
  theme_bw() +
  theme(legend.position = "bottom")  



## Box-plot with jitter points by sample 
ggplot(df, aes(x = Sample, y = LogSumGene, fill = Sample)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +  # Boxplot without outliers for clarity
  geom_jitter(shape = 21, size = 0.5, alpha = 0.5, width = 0.2) +  # Add jittered points
  scale_fill_brewer(palette = "Set2") +  # Nice color palette
  labs(title = "Log-Transformed Sum Gene Counts per Sample",
       x = "Sample",
       y = "Log10(Sum Gene Counts + 1)") +
  theme_bw() +
  theme(axis.text = element_text(size = 8, angle = 90)) +
  theme(legend.position = "none") 


dev.off()



# # Plot histogram by sample
# ggplot(df, aes(x = SumUMI, fill = Sample)) +
#   geom_histogram(binwidth = 500, alpha = 0.6, position = "identity") +
#   labs(title = "Distribution of Sum UMI Counts per Spot by Sample",
#        x = "Sum UMI Counts",
#        y = "Frequency") + 
#   theme_bw() + 
#   theme(legend.position="none") +
#   theme(axis.text = element_text(size = 6)) +
#   facet_wrap(~ Sample, scales = "free_y")  # Separate plots per sample
# 
# # Plot density distribution by sample in one plot
# ggplot(df, aes(x = SumUMI, color = Sample, fill = Sample)) +
#   geom_density(alpha = 0.3) +
#   labs(title = "Density Distribution of Sum UMI Counts per Spot by Sample",
#        x = "Sum UMI Counts",
#        y = "Density") +
#   theme_bw()



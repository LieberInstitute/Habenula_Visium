#   For a specific value of k (looped over using an array), gather FICTURE
#   results to ultimately form cell-level extracellular scores for each cluster,
#   and export a CSV of these scores. This intermediate step (prior to any
#   interesting analysis) is decently expensive in runtime and memory, hence
#   the dedication of a full separate script to this task

library(here)
library(tidyverse)
library(data.table)
library(spatialLIBD)
library(HDF5Array)
library(sessioninfo)

k = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

spe_dir = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'spe_norm_filtered'
)
ficture_path = here(
        'processed-data', '10_HD_bin_level', 'probe_fix', 'cell_environment',
        'ficture_outputs', 'normalized', 'k_%d', 'analysis', 'nF%d.d_12',
        'normalized_joined_input.tsv.gz'
    ) |>
    sprintf(k, k)
extra_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'cell_environment',
    'extracellular_bins.csv.gz'
)
out_path = here(
        'processed-data', '10_HD_bin_level', 'probe_fix', 'cell_environment',
        'cell_profiles', 'k_%d.csv.gz'
    ) |>
    sprintf(k)
plot_dir = here('plots', '10_HD_bin_level', 'probe_fix', 'cell_environment')
ficture_colnames = c(
    'sample_id', 'barcode', 'factor_K1', 'factor_K2', 'factor_K3', 'factor_P1',
    'factor_P2', 'factor_P3'
)

dir.create(dirname(out_path), showWarnings = FALSE)

################################################################################
#   Read in extracellular bins and FICTURE results
################################################################################

message(Sys.time(), ' | Reading in FICTURE clusters...')
ficture_df = fread(ficture_path, select = ficture_colnames) |>
    as_tibble()

message(Sys.time(), ' | Reading in extracellular bins...')
extra_df = fread(extra_path) |>
    as_tibble() |>
    dplyr::rename(barcode = bin_id)

message('Proportion of bins dropped by FICTURE (by sample):')
extra_df |>
    inner_join(ficture_df, by = c('sample_id', 'barcode'), multiple = 'any') |>
    group_by(sample_id) |>
    summarize(prop_missing = mean(is.na(factor_K1))) |>
    ungroup() |>
    print()

message(Sys.time(), ' | Joining and computing cell-level scores...')
extra_df = extra_df |>
    #   This is much faster than first taking unique combinations of sample_id
    #   and barcode before joining
    inner_join(ficture_df, by = c('sample_id', 'barcode'), multiple = 'any') |>
    filter(!is.na(factor_K1))

################################################################################
#   Compute cell-level scores for each FICTURE cluster
################################################################################

#   Form score columns for each cluster based on adding up posterior
#   probabilities for the top 3 factors
for (this_k in seq_len(k) - 1) {
    extra_df[[paste0('score_', this_k)]] = (
        (extra_df$factor_K1 == this_k) * extra_df$factor_P1 +
        (extra_df$factor_K2 == this_k) * extra_df$factor_P2 +
        (extra_df$factor_K3 == this_k) * extra_df$factor_P3
    )
}

#   Add up scores for each cell
extra_df = extra_df |>
    #   For each cell and sample, sum up scores
    group_by(sample_id, cell_id) |>
    summarize(
        across(matches('^score_'), sum),
        num_bins = n()
    ) |>
    ungroup() |>
    #   Normalize scores so they add to 1 across all clusters
    mutate(
        temp_sum = rowSums(across(matches('^score_'))),
        key = paste(cell_id, sample_id, sep = '_')
    ) |>
    mutate(across(matches('^score_'), function(x) x / temp_sum)) |>
    select(key, num_bins, matches('^score_'))

################################################################################
#   Join scores with SPE and export a minimal tibble
################################################################################

message(Sys.time(), ' | Loading cell-level SPE...')
spe = loadHDF5SummarizedExperiment(spe_dir)

extra_df = tibble(
        key = spe$key, sample_id = spe$sample_id,
        cell_category = spe$labels_joint_source
    ) |>
    left_join(extra_df, by = 'key')

message('Proportion of cells missing extracellular clustered bins (overall):')
extra_df |>
    group_by(cell_category) |>
    summarize(prop_missing = mean(is.na(score_0))) |>
    ungroup() |>
    print()

message('Proportion of cells missing extracellular clustered bins (by sample):')
extra_df |>
    group_by(sample_id, cell_category) |>
    summarize(prop_missing = mean(is.na(score_0))) |>
    ungroup() |>
    print(n = 10)

if (k == 4) {
    #   Just to verify cell IDs and joining are correct, plot scores for the
    #   white-matter cluster spatially on Br9090
    spe$score_3 = extra_df$score_3

    p = vis_gene(
        spe, sampleid = 'H1-W369TJK_D1_9090', geneid = 'score_3',
        is_stitched = TRUE, point_size = 10, spatial = TRUE
    )
    
    png(file.path(plot_dir, 'Br9090_WM_k4.png'), width = 1000, height = 1000)
    print(p)
    dev.off()

    #   Demonstrate the sample-specific effect induced by the poor-quality H&E
    #   images in the last 2 samples: "cells" are so dense that many have no
    #   extracellular bins
    p = extra_df |>
        group_by(sample_id, cell_category) |>
        summarize(prop_remain = mean(!is.na(score_0))) |>
        ungroup() |>
        mutate(
            sample_id = factor(
                str_replace(sample_id, '.*_', 'Br'),
                levels = c('Br9090', 'Br8433', 'Br8667', 'Br8518', 'Br9037')
            )
        ) |>
        ggplot(aes(x = sample_id, y = prop_remain, fill = sample_id)) +
            geom_col() +
            facet_wrap(~cell_category) +
            theme_bw(base_size = 20) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
            labs(x = 'Sample ID', y = 'Prop. Cells with Surrounding Bins') +
            guides(fill = 'none')
    pdf(file.path(plot_dir, 'missing_cells.pdf'))
    print(p)
    dev.off()
}

extra_df |>
    filter(!is.na(score_0)) |>
    select(key, num_bins, matches('^score_')) |>
    write_csv(out_path)

session_info()

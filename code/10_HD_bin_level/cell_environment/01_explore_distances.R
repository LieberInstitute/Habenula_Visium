library(here)
library(tidyverse)
library(segmented)
library(sessioninfo)

occupation_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'cell_environment',
    'occupation.csv'
)

occupation_df = read_csv(occupation_path, show_col_types = FALSE)
model = lm(occupation ~ expansion_distance, data = occupation_df)

message('Segmented model summary (occupation fraction vs. expansion distance):')
model |>
    segmented(
        seg.Z = ~expansion_distance, psi = list(expansion_distance = 6)
    ) |>
    summary()

session_info()

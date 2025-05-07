import scanpy as sc
import os
from pyhere import here
import session_info
import datetime

import extracellular_bins_functions as ebf

sample_id_path = here('raw-data', 'sample_info', 'hd_sample_list.txt')
with open(sample_id_path, 'r') as f:
    all_samples = f.read().splitlines()

# sample_id = all_samples[int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1]
sample_id = all_samples[0]

pre_out_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix',
    f'{sample_id}_pre_bin2cell.h5ad'
)
plot_dir = here('plots', '10_HD_bin_level', 'probe_fix', 'cell_environment')
mpp = 0.3
random_state = 0
num_random_cells = 5
min_bins_per_cell = 4

os.makedirs(plot_dir, exist_ok=True)

adata = sc.read(pre_out_path)
adata = ebf.find_microenvironment(adata, expansion_distance = 4)
adata = ebf.drop_bad_secondary_cells(adata, min_bins_per_cell = min_bins_per_cell)
extracellular_df = ebf.export_and_plot(adata, plot_dir, sample_id, mpp)

session_info.show()

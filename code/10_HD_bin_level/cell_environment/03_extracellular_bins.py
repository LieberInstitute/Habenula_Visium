import scanpy as sc
import os
from pyhere import here
import session_info
import datetime
import pandas as pd

import extracellular_bins_functions as ebf

sample_id_path = here('raw-data', 'sample_info', 'hd_sample_list.txt')
plot_dir = here(
    'plots', '10_HD_bin_level', 'probe_fix', 'cell_environment', 'random_cells'
)
out_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'cell_environment',
    'extracellular_bins.csv.gz'
)

mpp = 0.3
min_bins_per_cell = 4
expansion_distance = 6

os.makedirs(plot_dir, exist_ok=True)

with open(sample_id_path, 'r') as f:
    all_samples = f.read().splitlines()

extracellular_df_list = []
for sample_id in all_samples:
    print(f"{datetime.datetime.now()} | Processing sample {sample_id}")

    pre_out_path = here(
        'processed-data', '09_HD_cell_level', 'probe_fix',
        f'{sample_id}_pre_bin2cell.h5ad'
    )
    adata = sc.read(pre_out_path)

    adata = ebf.find_microenvironment(
        adata, expansion_distance = expansion_distance
    )
    adata = ebf.drop_bad_secondary_cells(
        adata, min_bins_per_cell = min_bins_per_cell
    )
    extracellular_df_list.append(
        ebf.export_and_plot(adata, plot_dir, sample_id, mpp)
    )

print(f"{datetime.datetime.now()} | Merging and exporting")
extracellular_df = pd.concat(extracellular_df_list, axis = 0)
extracellular_df.to_csv(out_path, index = False)

session_info.show()

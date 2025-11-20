import scanpy as sc
import os
from pyhere import here
import session_info
import datetime
import pandas as pd

import extracellular_bins_functions as ebf

sample_id_path = here('raw-data', 'sample_info', 'hd_basic_info.csv')
plot_dir = here(
    'plots', '10_HD_bin_level', 'new_samples', 'cell_environment', 'random_cells'
)
out_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples', 'cell_environment',
    'extracellular_bins.csv.gz'
)

mpp = 0.3
min_bins_per_cell = 4
expansion_distance = 5

os.makedirs(plot_dir, exist_ok=True)

sample_info = pd.read_csv(sample_id_path)
all_samples = sample_info['sample_id'].tolist()

extracellular_df_list = []
for sample_id in all_samples:
    print(f"{datetime.datetime.now()} | Processing sample {sample_id}")

    pre_out_path = here(
        'processed-data', '09_HD_cell_level', 'new_samples',
        f'{sample_id}_pre_bin2cell.h5ad'
    )
    adata = sc.read(pre_out_path)

    adata = ebf.find_microenvironment(
        adata, expansion_distance = expansion_distance
    )

    #   Add 'cell_component' column for informative coloring of plots
    adata.obs['cell_component'] = 'Unlabeled'
    adata.obs.loc[adata.obs['microenvironment_secondary'] != 0, 'cell_component'] = 'Sec. Extracellular'
    adata.obs.loc[adata.obs['labels_gex'] != 0, 'cell_component'] = 'Sec. Cell Body'
    adata.obs.loc[adata.obs['microenvironment_primary'] != 0, 'cell_component'] = 'Prim. Extracellular'
    adata.obs.loc[adata.obs['labels_he_expanded'] != 0, 'cell_component'] = 'Prim. Cell Body'
    adata.obs.loc[adata.obs['labels_he'] != 0, 'cell_component'] = 'Prim. Nucleus'

    #   Surprisingly, bin2cell passively drops many secondary cells by giving all of
    #   their bins (either this or none of the bins) a 'none' value in
    #   'labels_joint_source'. Just drop such bins now (cell body only)
    adata = adata[
        (adata.obs['labels_joint_source'] == 'secondary') |
        (adata.obs['cell_component'] != 'Sec. Cell Body')
    ]

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

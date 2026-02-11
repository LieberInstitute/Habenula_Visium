#   Export a CSV of extracellular bins and their associated cells; export a
#   single dataset-wide AnnData of "cells" formed by aggregating extracellular
#   expression

import scanpy as sc
import os
from pyhere import here
import session_info
import datetime
import pandas as pd
import bin2cell as b2c
import anndata as ad

import extracellular_bins_functions as ebf

sample_id_path = here('raw-data', 'sample_info', 'hd_basic_info.csv')
plot_dir = here(
    'plots', '10_HD_bin_level', 'no_secondary', 'cell_environment', 'random_cells'
)
bin_out_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'extracellular_bins.csv.gz'
)
adata_out_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'extracellular_adata.h5ad'
)

mpp = 0.3
expansion_distance = 7

os.makedirs(plot_dir, exist_ok=True)

sample_info = pd.read_csv(sample_id_path)
all_samples = sample_info['sample_id'].tolist()

extracellular_df_list = []
adata_list = []
for sample_id in all_samples:
    print(f"{datetime.datetime.now()} | Processing sample {sample_id}")

    pre_out_path = here(
        'processed-data', '09_HD_cell_level', 'new_samples2',
        f'{sample_id}_pre_bin2cell.h5ad'
    )
    adata = sc.read(pre_out_path)

    ebf.find_microenvironment(
        adata, expansion_distance = expansion_distance
    )

    #   Add 'cell_component' column for informative coloring of plots
    adata.obs['cell_component'] = 'Unlabeled'
    adata.obs.loc[adata.obs['microenvironment_primary'] != 0, 'cell_component'] = 'Prim. Extracellular'
    adata.obs.loc[adata.obs['labels_he_expanded'] != 0, 'cell_component'] = 'Prim. Cell Body'
    adata.obs.loc[adata.obs['labels_he'] != 0, 'cell_component'] = 'Prim. Nucleus'

    extracellular_df_list.append(
        ebf.export_and_plot(adata, plot_dir, sample_id, mpp)
    )

    #   Drop intracellular bins
    adata = adata[adata.obs['cell_component'] == 'Prim. Extracellular', :]

    adata = b2c.bin_to_cell(
        adata, labels_key="microenvironment_joint",
        spatial_keys=["spatial", "spatial_cropped_150_buffer"]
    )

    #   Add key which matches the ordinary cell-level anndata (not
    #   extracellular)
    adata.obs['key'] = adata.obs.index + '_' + adata.obs['sample_id']
    adata.obs.set_index('key', inplace=True)

    adata_list.append(adata)

print(f"{datetime.datetime.now()} | Merging and exporting")
extracellular_df = pd.concat(extracellular_df_list, axis = 0)
extracellular_df.to_csv(bin_out_path, index = False)

adata = ad.concat(adata_list, axis=0)
sc.write(adata_out_path, adata)

session_info.show()

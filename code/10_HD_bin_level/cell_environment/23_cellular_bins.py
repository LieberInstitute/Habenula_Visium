#   We'll want to plot FICTURE clusters and Banksy cell types in the same plots.
#   We already exported the set of extracellular bins (currently
#   03_extracellular_bins.py), and now this script does the cellular ones

import scanpy as sc
from pyhere import here
import session_info
import datetime
import pandas as pd

sample_id_path = here('raw-data', 'sample_info', 'hd_basic_info.csv')
out_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'cellular_bins.csv.gz'
)

sample_info = pd.read_csv(sample_id_path)
all_samples = sample_info['sample_id'].tolist()

cellular_df_list = []
for sample_id in all_samples:
    print(f"{datetime.datetime.now()} | Processing sample {sample_id}")

    pre_out_path = here(
        'processed-data', '09_HD_cell_level', 'new_samples2',
        f'{sample_id}_pre_bin2cell.h5ad'
    )
    adata = sc.read(pre_out_path)

    adata = adata[adata.obs['labels_he_expanded'] != 0, :]

    adata.obs['cell_key'] = adata.obs['labels_he_expanded'].astype(str) + '_' + sample_id
    adata.obs['bin_id'] = adata.obs.index + '_' + sample_id
    cellular_df_list.append(adata.obs[['bin_id', 'cell_key']])

print(f"{datetime.datetime.now()} | Merging and exporting")
cellular_df = pd.concat(cellular_df_list, axis = 0)
cellular_df.to_csv(out_path, index = False)

session_info.show()

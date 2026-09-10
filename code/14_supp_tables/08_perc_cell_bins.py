import scanpy as sc
import os
from pyhere import here
import session_info
import pandas as pd

sample_info_path = here('raw-data', 'sample_info', 'hd_basic_info.csv')

sample_info = pd.read_csv(sample_info_path)
all_samples = sample_info['sample_id'].tolist()

n_cellular = 0
n_total = 0

for sample_id in all_samples:
    pre_out_path = here(
        'processed-data', '09_HD_cell_level', 'new_samples2',
        f'{sample_id}_pre_bin2cell.h5ad'
    )
    adata = sc.read(pre_out_path)

    labels = adata.obs['labels_he_expanded']
    n_total += len(labels)
    n_cellular += (labels != 0).sum()

perc_cellular = 100 * n_cellular / n_total
print(f'{perc_cellular:.2f}% of bins are cellular ({n_cellular}/{n_total})')

session_info.show()

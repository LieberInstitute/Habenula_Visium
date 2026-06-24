#   Export an AnnData for each donor with one "cell" per primary nucleus

import os
from pyhere import here
import session_info
import bin2cell as b2c
import anndata as ad
import pandas as pd

sample_info_path = here('raw-data', 'sample_info', 'hd_basic_info.csv')
sample_info = pd.read_csv(sample_info_path)
task_id = int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1
sample_id = sample_info.iloc[task_id]['sample_id']

in_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2',
    f'{sample_id}_pre_bin2cell.h5ad'
)
out_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'contamination',
    'bin2cell', f'{sample_id}.h5ad'
)

os.makedirs(os.path.dirname(out_path), exist_ok=True)

adata = ad.read_h5ad(in_path)

#   Assumptions (which I interactively tested) that ensure the resulting cell IDs
#   match the cell-level object used in the main Visium HD analyses:
#   - In all bins where 'labels_he' in nonzero, 'labels_joint' shares the same label
#   - After bin_to_cell(), 'object_id' matches the label for primary cells (which I
#     tested by checking that bin_count always matches the number of rows for the
#     corresponding ID)
adata = b2c.bin_to_cell(
    adata, labels_key="labels_he",
    spatial_keys=["spatial", "spatial_cropped_150_buffer"]
)

adata.write_h5ad(out_path)

session_info.show()

import matplotlib.pyplot as plt
import scanpy as sc
import os
from pyhere import here
import session_info
import bin2cell as b2c
import datetime

prob_tag = "default"

sample_id_path = here('raw-data', 'sample_info', 'hd_sample_list.txt')
with open(sample_id_path, 'r') as f:
    all_samples = f.read().splitlines()
sample_id = all_samples[int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1]

stardist_dir = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'refine_segmentations',
    'stardist'
)
out_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'refine_segmentations',
    'adata', f'{sample_id}_{prob_tag}.h5ad'
)
orig_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix',
    f'{sample_id}_pre_bin2cell.h5ad'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'probe_fix', 'refine_segmentations', prob_tag
)
mpp = 0.3
num_random_cells = 5

os.makedirs(stardist_dir, exist_ok=True)
os.makedirs(plot_dir, exist_ok=True)
os.makedirs(out_path.parent, exist_ok=True)

################################################################################
#   Perform nuclear-based ("primary") segmentation
################################################################################

adata = sc.read(orig_path)

print(f"{datetime.datetime.now()} | Performing nuclear-based ('primary') segmentation")

#   Segment nuclei on H&E image
b2c.stardist(
    image_path=os.path.join(
        stardist_dir, f'{sample_id}_{prob_tag}.tiff'
    ),
    labels_npz_path=os.path.join(
        stardist_dir, f'{sample_id}_{prob_tag}.npz'
    ),
    stardist_model="2D_versatile_he"
)

#   Add segmentations to object
b2c.insert_labels(
    adata, 
    labels_npz_path=os.path.join(
        stardist_dir, f'{sample_id}_{prob_tag}.npz'
    ), 
    basis="spatial", 
    spatial_key="spatial_cropped_150_buffer",
    mpp=mpp, 
    labels_key="labels_primary"
)

sc.write(out_path, adata)

session_info.show()

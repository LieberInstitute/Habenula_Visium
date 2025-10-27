#   At the cell-level QC stage, it became clear that H1-6FX4YN3_A1_3942 and
#   H1-6FX4YN3_D1_9902 had striping effects in the positions of cells.
#   Currently, this effect is known in Visium HD, but the recommended correction
#   generally has worsened the striping effect (https://github.com/Teichlab/bin2cell/issues/45).
#   Here we'll test for these samples specifically if "correction" is overall
#   beneficial.

import matplotlib.pyplot as plt
import scanpy as sc
import os
from pyhere import here
import session_info
import bin2cell as b2c
import datetime

import pandas as pd
import numpy as np

sample_info_path = here('raw-data', 'sample_info', 'hd_basic_info.csv')
sample_info = pd.read_csv(sample_info_path)
task_id = int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1
sample_id = sample_info.iloc[task_id]['sample_id']
spaceranger_dir = sample_info.iloc[task_id]['spaceranger_dir']

stardist_dir = here(
    'processed-data', '09_HD_cell_level', 'new_samples', 'stardist'
)
temp_dir = '/fastscratch/myscratch/neagles'
sr_dir = here(
    spaceranger_dir, 'outs', 'binned_outputs', 'square_002um'
)
sr_spatial_dir = here(
    spaceranger_dir, 'outs', 'spatial'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'new_samples', 'bin2cell', 'stripe_test'
)
raw_image_path = here('raw-data', 'images', 'vis-hd', f'{sample_id}.tif')
mpp = 0.3

os.makedirs(stardist_dir, exist_ok=True)
os.makedirs(plot_dir, exist_ok=True)

print(f"{datetime.datetime.now()} | Building and preprocessing AnnData")

#   In the tutorial at https://nbviewer.org/github/Teichlab/bin2cell/blob/main/notebooks/demo.ipynb,
#   the filtered feature matrix, with additional gene-filtering steps, is used
#   before using specific settings to segment the gene-expression-based image
#   ('secondary segmentation'). While we're interested in retaining all genes
#   (i.e. using the raw feature matrix as in 'adata'), we want similar secondary
#   segmentation behavior as in the tutorial, hence 'adata_filtered'
adata_filtered = b2c.read_visium(
    sr_dir,
    count_file = 'filtered_feature_bc_matrix.h5',
    source_image_path = raw_image_path,
    spaceranger_image_path = sr_spatial_dir
)

#   Create a scaled H&E image attached to the object
b2c.scaled_he_image(
    adata_filtered,
    mpp = mpp,
    save_path = os.path.join(temp_dir, f'he_{sample_id}.tiff')
)

#   Require bins with nonzero counts
sc.pp.filter_cells(adata_filtered, min_counts=1)
sc.pp.filter_genes(adata_filtered, min_cells=3)

#   Attempt destriping, tracking UMI counts before and after
adata_filtered.obs['sum_umi_before'] = adata_filtered.X.sum(axis = 1)
b2c.destripe(adata_filtered)
adata_filtered.obs['sum_umi_after'] = adata_filtered.X.sum(axis = 1)

#   Cap UMI at the top for a more dynamic color scale in plots
adata_filtered.obs['capped_umi_before'] = np.clip(
    adata_filtered.obs['sum_umi_before'],
    0,
    (
        np.median(adata_filtered.obs['sum_umi_before']) + 
        np.std(adata_filtered.obs['sum_umi_before'])
    )
)
adata_filtered.obs['capped_umi_after'] = np.clip(
    adata_filtered.obs['sum_umi_after'],
    0,
    (
        np.median(adata_filtered.obs['sum_umi_after']) + 
        np.std(adata_filtered.obs['sum_umi_after'])
    )
)

#   Take a small subset in the top-right corner of the data
scale_factor = 0.3
mask = ( 
    (adata_filtered.obs['array_row'] >= (1 - scale_factor) * adata_filtered.obs['array_row'].max()) & 
    (adata_filtered.obs['array_col'] >= (1 - scale_factor) * adata_filtered.obs['array_col'].max())
)
adata_small = adata_filtered[mask, :]

#   Plot UMI counts before and after
for color_var in ['capped_umi_before', 'capped_umi_after']:
    sc.pl.spatial(
        adata_small, color=color_var, img_key=f"{mpp}_mpp_150_buffer",
        basis="spatial_cropped_150_buffer"
    )
    plt.savefig(
        os.path.join(plot_dir, f'{sample_id}_{color_var}.png')
    )
    plt.close('all')

session_info.show()

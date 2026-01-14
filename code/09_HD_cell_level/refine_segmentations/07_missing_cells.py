#   How often does something significantly bright show up in the gene-expression
#   image that is not captured by nuclear segmentations? In other words, is
#   secondary segmentation even needed (and likely to not just pick up noise)?

import matplotlib.pyplot as plt
import scanpy as sc
import os
from pyhere import here
import session_info
import bin2cell as b2c
import pandas as pd
import numpy as np
import scipy
from PIL import Image
import skimage

task_id = int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1

sample_info_path = here('raw-data', 'sample_info', 'hd_basic_info.csv')
sample_info = pd.read_csv(sample_info_path)
sample_id = sample_info.iloc[task_id]['sample_id']
spaceranger_dir = sample_info.iloc[task_id]['spaceranger_dir']

stardist_dir = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'stardist'
)
pre_out_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2',
    f'{sample_id}_pre_bin2cell.h5ad'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'new_samples2', 'refine_segmentations',
    'missing_cells'
)
mpp = 0.3
num_images = 2
image_width = 100

os.makedirs(plot_dir, exist_ok=True)

rng = np.random.default_rng(seed=0)

################################################################################
#   Functions
################################################################################

#   In this script we want to crop the labels and image separately (like putting
#   primary segmentations on top of the gene-expression image). Reimplement
#   b2c.view_labels for this purpose. Based on https://github.com/Teichlab/bin2cell/blob/9845c04091a3dea5db6b6b3f94f71431f4f554da/bin2cell/bin2cell.py#L248C5-L299C15
def custom_view_labels(image_path, labels_npz_path, crop_image, crop_labels, stardist_normalize=False):
    labels_sparse = scipy.sparse.load_npz(labels_npz_path)
    #determine memory efficient dtype to load the image as
    #if we'll be normalising, we want np.float16 for optimal RAM footprint
    #otherwise use np.uint8
    if stardist_normalize:
        dtype = np.float16
    else:
        dtype = np.uint8
    
    #PIL is better at handling crops memory efficiently than cv2
    img = Image.open(image_path)
    #ensure that it's in RGB (otherwise there's a single channel for greyscale)
    img = np.array(img.crop(crop_image).convert('RGB'), dtype=dtype)
    #subset labels to area of interest
    #crop is (left, upper, right, lower)
    #https://pillow.readthedocs.io/en/stable/reference/Image.html#PIL.Image.Image.crop
    #upper:lower, left:right
    labels_sparse = labels_sparse[crop_labels[1]:crop_labels[3], crop_labels[0]:crop_labels[2]]
    labels_sparse = labels_sparse.tocoo()
    
    #optionally normalise image
    if stardist_normalize:
        img = b2c.normalize(img)
        #actually cap the values - currently there are sub 0 and above 1 entries
        img[img<0] = 0
        img[img>1] = 1
        #turn back to uint8 for internal consistency
        img = (255*img).astype(np.uint8)
    
    border_sparse = scipy.sparse.coo_matrix(
        skimage.segmentation.find_boundaries(np.array(labels_sparse.todense()))
    )
    #   Not sure why the dimensions are slightly too big...
    assert border_sparse.row.max() <= img.shape[0] * 1.05
    assert border_sparse.col.max() <= img.shape[1] * 1.05
    border_sparse.row[border_sparse.row >= img.shape[0]] = img.shape[0] - 1
    border_sparse.col[border_sparse.col >= img.shape[1]] = img.shape[1] - 1

    #can now easily colour the borders similar to what was done for the fill
    img[border_sparse.row, border_sparse.col, :] = [255, 255, 0]
    
    return img

################################################################################
#   Plot primary segmentations over gene-expression image
################################################################################

adata = sc.read_h5ad(pre_out_path)

#   Sample random primary cells, which will just decide the centers of images
#   to plot
cell_indices = rng.choice(
    adata.obs['labels_he'].unique(), num_images, replace = False
)

for i, cell_index in enumerate(cell_indices):
    #   Determine region for plots
    array_row_center = adata.obs.loc[
        adata.obs['labels_he'] == cell_index, 'array_row'
    ].mean()
    array_col_center = adata.obs.loc[
        adata.obs['labels_he'] == cell_index, 'array_col'
    ].mean()

    mask = (
        (adata.obs['array_row'] >= array_row_center - image_width // 2) & 
        (adata.obs['array_row'] <= array_row_center + image_width // 2) & 
        (adata.obs['array_col'] >= array_col_center - image_width // 2) & 
        (adata.obs['array_col'] <= array_col_center + image_width // 2)
    )

    print(f'Image {i+1} bounds: array_row:({array_row_center - image_width // 2}, {array_row_center + image_width // 2}), array_col:({array_col_center - image_width // 2}, {array_col_center + image_width // 2})')

    #   Plot primary segmentations over gene-expression image
    crop_primary = b2c.get_crop(
        adata[mask], basis="spatial", spatial_key="spatial_cropped_150_buffer",
        mpp=mpp
    )
    crop_secondary = b2c.get_crop(adata[mask], basis="array", mpp=mpp)
    rendered = custom_view_labels(
        image_path = os.path.join(stardist_dir, f'gex_{sample_id}.tiff'),
        labels_npz_path = os.path.join(stardist_dir, f'he_{sample_id}.npz'),  
        crop_image = crop_secondary, crop_labels = crop_primary,
        stardist_normalize = True
    )
    plt.imshow(rendered)
    plt.savefig(os.path.join(plot_dir, f'{sample_id}_{i+1}_secondary.png'))
    plt.close('all')

session_info.show()

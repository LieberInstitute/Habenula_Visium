#   Is there a difference in secondary segmentation quality between
#   ambiguous and non-ambiguous Banksy clusters? Plot around a random
#   secondary cell from clusters 4 and 27, as well as one cell from any other
#   cluster, for each sample.

import matplotlib.pyplot as plt
import scanpy as sc
import os
from pyhere import here
import session_info
import bin2cell as b2c
import pandas as pd

sample_info_path = here('raw-data', 'sample_info', 'hd_basic_info.csv')
sample_info = pd.read_csv(sample_info_path)

ad_paths = {
    sample_id: here(
        'processed-data', '09_HD_cell_level', 'new_samples2',
        f'{sample_id}.h5ad'
    )
    for sample_id in sample_info['sample_id']
}
pre_ad_paths = {
    sample_id: here(
        'processed-data', '09_HD_cell_level', 'new_samples2',
        f'{sample_id}_pre_bin2cell.h5ad'
    )
    for sample_id in sample_info['sample_id']
}
banksy_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'banksy', 'lambda0_2',
    'leiden_res1_7.csv'
)
stardist_dir = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'stardist'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'new_samples2', 'secondary_segmentations'
)
mpp = 0.3

os.makedirs(plot_dir, exist_ok=True)

banksy_df = pd.read_csv(banksy_path, index_col='key')

#   Collect one image per sample ID per cluster (ambiguous and non-ambiguous)
rendered_images = []
for sample_id in sample_info['sample_id']:
    #   Load in just the secondary cells
    ad = sc.read_h5ad(ad_paths[sample_id])
    ad = ad[ad.obs['labels_joint_source'] == 'secondary', :]
    
    ad_pre = sc.read_h5ad(pre_ad_paths[sample_id])
    
    #   Join in clustering results
    ad.obs.index = ad.obs.index + '_' + sample_id
    ad.obs['banksy'] = banksy_df['banksy_lambda0_2']
    ad = ad[~ad.obs['banksy'].isna(), :].copy()
    ad.obs['banksy'] = ad.obs['banksy'].astype(int)
    
    #   Create views for ambiguous clusters and a random non-ambiguous cluster
    ad_views = [
        ad[ad.obs['banksy'] == 4, :],
        ad[ad.obs['banksy'] == 27, :],
        ad[~ad.obs['banksy'].isin([4, 27]), :]
    ]
    
    sample_images = []
    for ad_view in ad_views:
        col_data = ad_view.obs.sample(n=1, random_state=0)
        
        mask = (
            (ad_pre.obs['array_row'] >= col_data['array_row'].values[0] - 25) & 
            (ad_pre.obs['array_row'] <= col_data['array_row'].values[0] + 25) & 
            (ad_pre.obs['array_col'] >= col_data['array_col'].values[0] - 25) & 
            (ad_pre.obs['array_col'] <= col_data['array_col'].values[0] + 25)
        )
        crop = b2c.get_crop(ad_pre[mask], basis="array", mpp=mpp)
        
        rendered = b2c.view_labels(
            image_path = os.path.join(
                stardist_dir, f'gex_{sample_id}.tiff'
            ),
            labels_npz_path = os.path.join(
                stardist_dir, f'gex_{sample_id}.npz'
            ),  
            crop = crop,
            stardist_normalize = True
        )
        sample_images.append(rendered)
    
    rendered_images.append(sample_images)

# Create grid of images
num_samples = len(rendered_images)
num_views = 3
fig, axes = plt.subplots(
    num_samples, num_views, figsize=(num_views * 5, num_samples * 5)
)

for i, (sample_id, sample_images) in enumerate(zip(sample_info['sample_id'].tolist(), rendered_images)):
    for j, rendered in enumerate(sample_images):
        axes[i, j].imshow(rendered)
        axes[i, j].axis('off')
        if i == 0:
            view_labels = ['Cluster 4', 'Cluster 27', 'Other Clusters']
            axes[i, j].set_title(view_labels[j])
        if j == 0:
            axes[i, j].set_ylabel(sample_id, rotation=90, size='large')

plt.tight_layout()
plt.savefig(
    os.path.join(plot_dir, 'secondary_segmentations_grid.png'),
    dpi=150, bbox_inches='tight'
)
plt.close('all')

session_info.show()

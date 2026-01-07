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
banksy_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'banksy', 'lambda0_2',
    'leiden_res1_7.csv'
)
stardist_dir = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'stardist'
)

banksy_df = pd.read_csv(banksy_path, index_col='key')

for sample_id in sample_info['sample_id']:
    #   Load in just the secondary cells
    ad = sc.read_h5ad(ad_paths[sample_id])
    ad = ad[ad.obs['labels_joint_source'] == 'secondary', :]

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
    for ad_view in ad_views:
        col_data = ad_view.obs.sample(n=1, random_state=0)

        mask = (
            (ad.obs['array_row'] >= col_data['array_row'].values[0] - 25) & 
            (ad.obs['array_row'] <= col_data['array_row'].values[0] + 25) & 
            (ad.obs['array_col'] >= col_data['array_col'].values[0] - 25) & 
            (ad.obs['array_col'] <= col_data['array_col'].values[0] + 25)
        )

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
    plt.imshow(rendered)
    plt.savefig(
        os.path.join(plot_dir, f'{sample_id}_secondary_segmentation{i+1}.png')
    )
    plt.close('all')

    
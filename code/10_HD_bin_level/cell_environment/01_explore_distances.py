#   We want to look at the microenvironment around cells, but what radius from
#   each cell centroid is appropriate? This script tries to find an optimum by
#   finding the point where the number of bins occupied by the smallest 50% of cells
#   stops increasing despite further expansion (indicating cells are too dense to
#   benefit much from further expansion).

import scanpy as sc
import os
from pyhere import here
import session_info
import matplotlib.pyplot as plt
import pandas as pd

import extracellular_bins_functions as ebf

sample_info_path = here('raw-data', 'sample_info', 'hd_basic_info.csv')
sample_info = pd.read_csv(sample_info_path)
sample_id = sample_info['sample_id'].iloc[
    int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1
]

ficture_cols = ['sample_id', 'barcode', 'FICTURE_k4']
ficture_WM_cluster = 3
mpp = 0.3

ficture_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
    'bin_level_clusters_batch.csv.gz'
)
pre_out_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2',
    f'{sample_id}_pre_bin2cell.h5ad'
)
df_out_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'occupation', f'{sample_id}.csv'
)
plot_dir = here('plots', '10_HD_bin_level', 'no_secondary', 'cell_environment')

os.makedirs(plot_dir, exist_ok=True)
os.makedirs(df_out_path.parent, exist_ok=True)

################################################################################
#   Determine optimal expansion distance
################################################################################

adata = sc.read(pre_out_path)

#   Number of bins occupied by each cell before expansion
num_cells = (
    adata.obs
        .loc[
            adata.obs['labels_joint'] != 0,
            'labels_joint'
        ]
        .value_counts()
        .sort_values()
)

expansion_distance = [0]
occupation = [num_cells[:int(len(num_cells) / 2)].sum() / adata.shape[0]]

for exp_d in list(range(1, 11)) + [15, 20]:
    adata = ebf.find_microenvironment(adata, expansion_distance = exp_d)
    #   Take the combined number of bins occupied by cells in the lower 50% by
    #   size. The idea is that expansion should not increase this number past a
    #   certain point due to high density of cells
    num_cells = (
        adata.obs
            .loc[
                adata.obs['microenvironment_primary'] != 0,
                'microenvironment_primary'
            ]
            .value_counts()
            .sort_values()
    )
    #
    expansion_distance.append(exp_d)
    occupation.append(
        num_cells[:int(len(num_cells) / 2)].sum() / adata.shape[0]
    )

occupation_df = pd.DataFrame(
    {
        'expansion_distance': expansion_distance,
        'occupation': occupation
    }
)
occupation_df.to_csv(df_out_path, index = False)

#   For this particular sample, there are white matter tracts where we'll
#   explore cell density
if sample_id == 'H1-W369TJK_D1_9090':
    #   Add in FICTURE clusters (k = 4) for this sample to adata.obs
    ficture_df = pd.read_csv(ficture_path, usecols = ficture_cols)
    ficture_df = (
        ficture_df
            .loc[
                ficture_df['sample_id'] == sample_id,
                ['barcode', 'FICTURE_k4']
            ]
            .set_index('barcode')
    )
    adata.obs['FICTURE_k4'] = ficture_df['FICTURE_k4']

    adata = ebf.find_microenvironment(adata, expansion_distance = 6)

    #   Label each bin with a cellular component
    adata.obs['cell_component'] = 'Unlabeled'
    adata.obs.loc[adata.obs['microenvironment_primary'] != 0, 'cell_component'] = 'Prim. Extracellular'
    adata.obs.loc[adata.obs['labels_he_expanded'] != 0, 'cell_component'] = 'Prim. Cell Body'
    adata.obs.loc[adata.obs['labels_he'] != 0, 'cell_component'] = 'Prim. Nucleus'

    ############################################################################
    #   Plot cells in a region rich in white matter
    ############################################################################

    #   For this sample, there are two slices of tissue that center at different
    #   values of array_col. Find the center of the white matter in one of the
    #   tissue slices
    small_obs = adata.obs[
        (adata.obs['array_col'] > adata.obs['array_col'].max() / 2) &
        (adata.obs['FICTURE_k4'] == ficture_WM_cluster)
    ]
    row_centroid = small_obs['array_row'].mean()
    col_centroid = small_obs['array_col'].mean()

    #   After trial and error, a bit of a pertubation from the centroid yields a
    #   region rich in white matter
    small_adata = adata[
        (adata.obs['array_row'] >= row_centroid - 100 - 40) &
        (adata.obs['array_row'] <= row_centroid - 100 + 40) &
        (adata.obs['array_col'] >= col_centroid - 40) &
        (adata.obs['array_col'] <= col_centroid + 40),
        :
    ]
    print(f'Distribution of FICTURE clusters in white-matter-rich region (WM is {ficture_WM_cluster}):')
    print(small_adata.obs['FICTURE_k4'].value_counts())

    small_adata.obs['FICTURE_k4'] = small_adata.obs['FICTURE_k4'].astype(str)
    sc.pl.spatial(
        small_adata, color=[None, "cell_component", "FICTURE_k4"],
        img_key=f"{mpp}_mpp_150_buffer", basis="spatial_cropped_150_buffer"
    )
    plt.savefig(
        os.path.join(plot_dir, f'{sample_id}_WM.png')
    )
    plt.close('all')

    ############################################################################
    #   Plot cells in a region without much white matter
    ############################################################################

    row_centroid = int(adata.obs['array_row'].max() / 2)
    col_centroid = int(3 * adata.obs['array_col'].max() / 4)

    small_adata = adata[
        (adata.obs['array_row'] >= row_centroid - 40) &
        (adata.obs['array_row'] <= row_centroid + 40) &
        (adata.obs['array_col'] >= col_centroid - 40) &
        (adata.obs['array_col'] <= col_centroid + 40),
        :
    ]
    print(f'Distribution of FICTURE clusters in white-matter-absent region (WM is {ficture_WM_cluster}):')
    print(small_adata.obs['FICTURE_k4'].value_counts())
    print('Distribution of full tissue sample:')
    print(adata.obs['FICTURE_k4'].value_counts())

    small_adata.obs['FICTURE_k4'] = small_adata.obs['FICTURE_k4'].astype(str)
    sc.pl.spatial(
        small_adata, color=[None, "cell_component", "FICTURE_k4"],
        img_key=f"{mpp}_mpp_150_buffer", basis="spatial_cropped_150_buffer"
    )
    plt.savefig(
        os.path.join(plot_dir, f'{sample_id}_nonWM.png')
    )
    plt.close('all')

session_info.show()

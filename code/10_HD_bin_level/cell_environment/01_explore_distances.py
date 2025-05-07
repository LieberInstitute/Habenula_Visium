#   We want to look at the microenvironment around cells, but what radius from
#   each cell centroid is appropriate? This script tries to find an optimum by
#   finding the point where the number of bins occupied by the smallest 50% of cells
#   stops increasing despite further expansion (indicating cells are too dense to
#   benefit much from further expansion).

import scanpy as sc
import os
from pyhere import here
import session_info
import datetime
import matplotlib.pyplot as plt
import plotnine as pn
import pandas as pd

import extracellular_bins_functions as ebf

sample_id_path = here('raw-data', 'sample_info', 'hd_sample_list.txt')
with open(sample_id_path, 'r') as f:
    all_samples = f.read().splitlines()

# sample_id = all_samples[int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1]
sample_id = all_samples[0]

pre_out_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix',
    f'{sample_id}_pre_bin2cell.h5ad'
)
plot_dir = here('plots', '10_HD_bin_level', 'probe_fix', 'cell_environment')

os.makedirs(plot_dir, exist_ok=True)

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

for exp_d in range(1, 11):
    adata = ebf.find_microenvironment(adata, expansion_distance = exp_d)
    #   Take the combined number of bins occupied by cells in the lower 50% by
    #   size. The idea is that expansion should not increase this number past a
    #   certain point due to high density of cells
    num_cells = (
        adata.obs
            .loc[
                adata.obs['microenvironment_joint'] != 0,
                'microenvironment_joint'
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

#   Plot fraction of bins occupied by the smallest 50% of cells against
#   expansion distance
(
    pn.ggplot(
            occupation_df, pn.aes(x = 'expansion_distance', y = 'occupation')
        ) +
        pn.geom_line() +
        pn.theme_bw(base_size = 20) +
        pn.labs(
            x = 'Expansion Distance (Num. Bins)',
            y = 'Fraction of Bins Occupied\nby Smallest 50%',
        )
)
plt.savefig(os.path.join(plot_dir, 'occupation.png'))
plt.close('all')

session_info.show()

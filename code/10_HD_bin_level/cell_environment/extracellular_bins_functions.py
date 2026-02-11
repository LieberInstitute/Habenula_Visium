#   Functions for 01_extracellular_bins.py

import bin2cell as b2c
import matplotlib.pyplot as plt
import scanpy as sc
import os
import pandas as pd

#   Find microenvironment around primary segmentations, expanding around the
#   cell bodies by [expansion_distance] bins. Modify the AnnData in place
def find_microenvironment(adata, expansion_distance):
    #   Label microenvironment around primary segmentations (nuclei)
    b2c.expand_labels(
        adata, 
        labels_key='labels_he_expanded', 
        expanded_labels_key="microenvironment_primary",
        max_bin_distance = expansion_distance
    )

    #   Sanity checks: cell labels should be preserved when expanding
    #   microenvironment. Microenvironment and cell labels should be one to one
    mask = adata.obs['labels_he_expanded'] != 0
    assert all(adata.obs['labels_he_expanded'][mask] == adata.obs['microenvironment_primary'][mask]), "Primary microenvironment labels don't always match their original cells"
    assert len(adata.obs['labels_he_expanded'].unique()) == len(adata.obs['microenvironment_primary'].unique()), "Primary microenvironment labels should be one to one with cell labels"

#   Visualize cell segmentations and surrounding microenvironment. Return a 
#   DataFrame containing just extracellular bins, ready for export
def export_and_plot(adata, plot_dir, sample_id, mpp, random_state = 0):
    #---------------------------------------------------------------------------
    #   Form DataFrame of extracellular bins for export
    #---------------------------------------------------------------------------

    extracellular_df = adata.obs[
        adata.obs['cell_component'].isin(
            ['Prim. Extracellular', 'Sec. Extracellular']
        )
    ].copy()

    #   For each secondary cell, find the corresponding label in 'labels_joint'
    label_key = (
        adata.obs
            .loc[adata.obs['labels_joint_source'] == 'secondary']
            .drop_duplicates(subset = 'labels_gex', keep = 'first')
            [['labels_gex', 'labels_joint']]
            .rename(
                {
                    'labels_joint': 'corresponding_joint',
                    'labels_gex': 'microenvironment_secondary'
                },
                axis = 1
            )
    )

    #   Add 'cell_id' column
    temp = extracellular_df.index
    extracellular_df = pd.merge(
        extracellular_df, label_key, how = 'left',
        on = 'microenvironment_secondary'
    )
    extracellular_df.index = temp
    extracellular_df['cell_id'] = extracellular_df['corresponding_joint']
    mask = extracellular_df['cell_component'] == 'Prim. Extracellular'
    extracellular_df.loc[mask, 'cell_id'] = extracellular_df.loc[
        mask, 'microenvironment_primary'
    ]
    assert all(~extracellular_df['cell_id'].isna())

    #   Clean up
    extracellular_df = (
        extracellular_df[['cell_id']]
            .reset_index(names = 'bin_id')
            .assign(sample_id = sample_id)
    )
    extracellular_df['cell_id'] = extracellular_df['cell_id'].astype(int)

    #---------------------------------------------------------------------------
    #   Visualize cell segmentations and surrounding microenvironment
    #---------------------------------------------------------------------------

    random_cell = (
        extracellular_df
            .loc[mask.values, :]
            .sample(n = 1, random_state = random_state)
            ['cell_id']
            .values[0]
    )
    small_adata = adata[adata.obs['microenvironment_primary'] == random_cell, :]

    small_adata = adata[
        (adata.obs['array_row'] >= small_adata.obs['array_row'].min() - 40) &
        (adata.obs['array_row'] <= small_adata.obs['array_row'].max() + 40) &
        (adata.obs['array_col'] >= small_adata.obs['array_col'].min() - 40) &
        (adata.obs['array_col'] <= small_adata.obs['array_col'].max() + 40),
        :
    ]

    sc.pl.spatial(
        small_adata, color=[None, "cell_component"],
        img_key=f"{mpp}_mpp_150_buffer", basis="spatial_cropped_150_buffer"
    )
    plt.savefig(
        os.path.join(plot_dir, f'{sample_id}_random_cells.png')
    )
    plt.close('all')

    return extracellular_df

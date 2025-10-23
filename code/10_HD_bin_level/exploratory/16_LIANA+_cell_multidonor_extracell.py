# Create the LIANA files for the multidonor analysis
import scanpy as sc
import decoupler as dc
import plotnine as p9
import liana as li
import numpy as np
import pandas as pd
import os
from pyhere import here
import matplotlib.pyplot as plt
import json
import numpy as np
from pyhere import here
import re
import session_info
from plotnine import *

#   Read input files
in_dir = here(
    'processed-data', '10_HD_bin_level', 'LIANA'
)
plot_dir= here(
    'processed-data', '10_HD_bin_level', 'LIANA', 'figure','habenula'
)
os.makedirs(plot_dir, exist_ok=True)

in_files = [ os.path.join(in_dir, f) for f in os.listdir(in_dir) if re.compile(r'.*\.h5ad$').match(f) ]
in_files = [f for f in in_files if "extracellular" in f]

#   Read in DataFrames of ligand-receptor stats for each donor and concatenate
lr_df_list = []
for f in in_files:
    donor_id = f.split('/')[-1].replace('.h5ad', '').replace('extracellular_lrdata_', '')
    lr_df = sc.read(f).var
    lr_df['donor_id'] = donor_id
    lr_df_list.append(lr_df)

lr_df = pd.concat(lr_df_list, axis=0, ignore_index=True)

min_num_donors=3
#   Require a pair to be present in some minimum number of donors
lr_df = lr_df[
    lr_df.groupby(['ligand', 'receptor'], observed=True)['ligand'].transform('count') >= min_num_donors
]

#   Average stats across samples
lr_mean_df = (
    lr_df
        .groupby(['ligand', 'receptor'], as_index=False)
        .agg({'mean': 'mean', 'morans': 'mean'})
)

lr_top_df = pd.concat(
    [
        lr_mean_df.sort_values("morans", ascending=False).head(),
        lr_mean_df.sort_values("mean", ascending=False).head()
    ]
)
lr_top_df = lr_top_df.drop_duplicates(subset=['ligand', 'receptor'])

top_pairs = list(lr_top_df['ligand'] + '^' + lr_top_df['receptor'])

#  first sample
for f in in_files:                      
    ad_lr = sc.read(f)                  
    sample_id = os.path.basename(f).replace(".h5ad", "").replace("extracellular_lrdata_", "")

    assert all([x in ad_lr.var.index for x in top_pairs])

    #   Fix spatial coordinates format
    ad_lr.obsm['spatial'] = np.asarray(ad_lr.obsm['spatial'])

    #   Plot scores, permutation-based p-values, and local categories
    fig = sc.pl.spatial(
        ad_lr, color=top_pairs, cmap='viridis', ncols=2,
        spot_size=80, show=False, return_fig=True
    )
    fig.savefig(os.path.join(plot_dir, f"extracellular_bivariate_scores_{sample_id}.pdf"), bbox_inches='tight')
    plt.close(fig)

    sc.pl.spatial(
        ad_lr, color=top_pairs, cmap='viridis_r', ncols=2,
        spot_size=80, layer='pvals'
    )
    plt.savefig(os.path.join(plot_dir, f"extracellular_bivariate_p_vals_{sample_id}.pdf"))
    plt.close('all')

    sc.pl.spatial(
        ad_lr, color=top_pairs, cmap='coolwarm', ncols=2,
        spot_size=80, layer='cats'
    )
    plt.savefig(os.path.join(plot_dir, f"extracellular_bivariate_cats_{sample_id}.pdf"))
    plt.close('all')

session_info.show()


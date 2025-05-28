import os
from pyhere import here
import session_info
import pandas as pd

k = int(os.getenv('SLURM_ARRAY_TASK_ID'))

extra_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'cell_environment',
    'cell_profiles', f'k_{k}.csv.gz'
)
banksy_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'banksy', 'lambda0_2',
    'leiden_res1_4.csv'
)
good_samples = [
    'H1-W369TJK_D1_9090', 'H1-MVPY9BW_A1_8433', 'H1-MVPY9BW_D1_8667'
]

extra_df = pd.read_csv(extra_path, index_col = 'key')
banksy_df = pd.read_csv(banksy_path, index_col = 'key')

#   Grab only good samples
mask = ['_'.join(x.split('_')[1:]) in good_samples for x in banksy_df.index]
banksy_df = banksy_df.loc[mask, :]
mask = ['_'.join(x.split('_')[1:]) in good_samples for x in extra_df.index]
extra_df = extra_df.loc[mask, :]

#   Add Banksy clusters to extra_df. All cells should have cluster assignments
extra_df['banksy'] = banksy_df['banksy_lambda0_2']
assert extra_df.isna().any().sum() == 0

#   Can we determine whether a cell is in the habenula or not based on its
#   extracellular FICTURE cluster distribution?

import os
from pyhere import here
import session_info
import pandas as pd
import numpy as np
from sklearn import tree, svm
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import classification_report
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
import matplotlib.pyplot as plt

k = int(os.getenv('SLURM_ARRAY_TASK_ID'))

extra_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'cell_environment',
    'cell_profiles', f'k_{k}.csv.gz'
)
banksy_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'banksy', 'lambda0_2',
    'leiden_res1_4.csv'
)
plot_path = here(
    'plots', '10_HD_bin_level', 'probe_fix', 'cell_environment', 'hb_tree.pdf'
)
good_samples = [
    'H1-W369TJK_D1_9090', 'H1-MVPY9BW_A1_8433', 'H1-MVPY9BW_D1_8667'
]
test_prop = 0.2
hb_clusters = [2, 11]
random_seed = 0
downsample_size = 5000

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

extra_df['is_hb'] = pd.Series(
    [
        'habenula' if x else 'not_habenula'
        for x in extra_df['banksy'].isin(hb_clusters)
    ],
    dtype = 'category', index = extra_df.index
)
extra_df['sample_id'] = pd.Series(
    ['_'.join(x.split('_')[1:]) for x in extra_df.index], dtype = 'category',
    index = extra_df.index
)

x_train, x_test, y_train, y_test = train_test_split(
    extra_df.filter(regex='^score_', axis = 1), extra_df['is_hb'],
    test_size = test_prop, random_state = random_seed,
    stratify = (extra_df['sample_id'].astype(str) + extra_df['is_hb'].astype(str)).astype('category')
)

#   Also downsample to make models train in a reasonable time
small_extra_df = extra_df.sample(
    n = int(downsample_size / (1 - test_prop)), random_state = random_seed
)
x_train_small, x_test_small, y_train_small, y_test_small = train_test_split(
    small_extra_df.filter(regex='^score_', axis = 1), small_extra_df['is_hb'],
    test_size = test_prop, random_state = random_seed,
    stratify = (small_extra_df['sample_id'].astype(str) + small_extra_df['is_hb'].astype(str)).astype('category')
)

################################################################################
#   Decision tree
################################################################################

#   A single decision tree, to see if we can get straightforward-to-interpret
#   results (rather than using an optimal model for classification)

model = tree.DecisionTreeClassifier(
    max_depth = 2 * k, min_samples_leaf = 0.05, random_state = random_seed,
    ccp_alpha = 0.001
)
model.fit(x_train, y_train)

print('---- Trying decision tree...')
print('Training report:\n', classification_report(y_train, model.predict(x_train)))
print('Test report:\n', classification_report(y_test, model.predict(x_test)))

plt.figure(figsize=(10, 5))
tree.plot_tree(
    model, class_names = model.classes_, feature_names = x_train.columns,
    filled = True, rounded = True, fontsize = 5
)
plt.savefig(plot_path)
plt.close('all')

################################################################################
#   SVM
################################################################################

#   Here we're trying to use a more powerful model to see what the best F1-score
#   we can get is (the idea being that the model performance puts an upper bound
#   on how well we can even predict habenula vs. not habenula from extracellular
#   environment)
tuned_parameters = [
    {
        'svc__kernel': ['rbf', 'poly'],
        'svc__C': np.logspace(-2, 2, 5)
    }
]

pipe = make_pipeline(
    StandardScaler(),
    svm.SVC(class_weight = 'balanced', random_state = random_seed)
)

#   Outputs should be integers
y_train_small_int = (y_train_small == 'habenula').astype('int')
y_test_small_int = (y_test_small == 'habenula').astype('int')

grid = GridSearchCV(pipe, tuned_parameters, cv = 5, scoring = 'f1')
grid.fit(x_train_small, y_train_small_int)

print('---- Trying SVM...')
print(
    'Training report:\n',
    classification_report(
        y_train_small_int, grid.best_estimator_.predict(x_train_small)
    )
)
print(
    'Test report:\n',
    classification_report(
        y_test_small_int, grid.best_estimator_.predict(x_test_small)
    )
)

session_info.show()

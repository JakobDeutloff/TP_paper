"""
Script to read data for the comparison of LH and random sampling (variance analysis) and the stability analysis of
our results
"""

import pandas as pd
import glob
import os
import numpy as np

# %% Load var analysis
path = r'Data/model_output/585_var_analysis/'
all_files_LH = glob.glob(os.path.join(path, r"*LH.csv"))
all_files_rand = glob.glob(os.path.join(path, r"*Rand.csv"))
li = []
for filename in all_files_LH:
    df = pd.read_csv(filename, index_col=[0], header=[0, 1, 2], parse_dates=[0])
    li.append(df)
keys = np.arange(len(all_files_LH))
keys = [str(x) for x in keys]
LH585 = pd.concat(li, axis=1, keys=keys)

li = []
for filename in all_files_rand:
    df = pd.read_csv(filename, index_col=[0], header=[0, 1, 2], parse_dates=[0])
    li.append(df)
keys = np.arange(len(all_files_LH))
keys = [str(x) for x in keys]
Rand585 = pd.concat(li, axis=1, keys=keys)

# %% Load stability analysis
path = r'Data/model_output/585_stab_analysis'
T = pd.read_csv(path + '/T.csv', index_col=0, parse_dates=[0], header=[0, 1])
T.columns.names = ['ens_nember', 'percentile']

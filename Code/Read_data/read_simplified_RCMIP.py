"""
Script to read the simplified RCMIP data used to run FaIR
"""

import pandas as pd

#%%
SSP_emms = pd.read_csv('Data/RCMIP_simplified/SSP_emms.csv', index_col=[0], header=[0, 1])
SSP_forc = pd.read_csv('Data/RCMIP_simplified/SSP_forc.csv', index_col=[0], header=[0, 1])

"""
This Script is used to calculate theyear where vertain percentiles of tipping probabilit are exceeded,
"""
# %%
from Code.Calibrate.find_distributions import P_perf_dists
from Code.Read_data.Read_SSP_output import read_SSP_outout, read_full_T
import pandas as pd
import numpy as np
import pickle

# %% Load data

#  Set ens name
ens_name = 'uncoupled_ensemble'
SSP_out = read_SSP_outout(ens_name)
T = read_full_T(ens_name)

# %% Build tipping probs dataframes

elements = ['PFTP', 'AMOC', 'GRIS', 'WAIS', 'AMAZ', 'BORF', 'AWSI', 'EAIS', 'EASB', 'GLCR', 'LABC', 'TUND', 'PFAT',
            'BARI', 'REEF', 'SAHL']
scenarios = SSP_out['T'].columns.levels[0]
percentile = SSP_out['T'].columns.levels[1]


# DF for years of crossing percentiles of total tipping probabilities
multicol = pd.MultiIndex.from_product([scenarios, elements])
tip_years = pd.DataFrame(index=percentile, columns=multicol)
tip_years.columns.names = ['scenario', 'element']

# DF for cumulative TEs from total probabilities
multicol = pd.MultiIndex.from_product([scenarios, percentile])
cum_TEs = pd.DataFrame(index=SSP_out['T'].index, columns=multicol, data=np.zeros([len(SSP_out['T'].index),
                                                                                  len(scenarios)*len(percentile)]))
cum_TEs.columns.names = ['scenario', 'percentile']

# DF for cum C emissions in GtC
tip_emms = pd.DataFrame(index=percentile, columns=multicol)
tip_emms.columns.names = ['scenario', 'element']

# DF for slow tipping probabilities
multicol = pd.MultiIndex.from_product([scenarios, elements])
slow_tip_probs = pd.DataFrame(index=[2450], columns=multicol)

# %% calculate year of exceedance for given percentiles
for ssp in tip_years.columns.levels[0].values:
    for perc in percentile:
        for element in elements:

            try:
                idx = np.where(SSP_out['Tip_probs_total'][ssp, element] > float(perc))[0][0]
                tip_years.loc[perc, (ssp, element)] = SSP_out['Tip_probs_total'].index[idx]
            except:
                pass

# %% Calculate cumulative C emissions at tipping years
for ssp in tip_years.columns.levels[0].values:
    # Cumulative carbon emissions (CO2 and CH4) in GtC
    cum_C = SSP_out['Emm'].loc[:, (ssp, 'carbon_dioxide')].cumsum() + SSP_out['Emm'].loc[:,
                                                                      (ssp, 'methane')].cumsum() / 1e3
    for perc in percentile:
        for element in elements:
            try:
                tip_emms.loc[perc, (ssp, element)] = cum_C[tip_years.loc[perc, (ssp, element)]]
            except:
                pass

# %% Calculate cumulative TEs that have crossed the probability of exceeding a given percentile
for ssp in tip_years.columns.levels[0].values:
    for perc in percentile:
        for elem in elements:
            year = tip_years.loc[perc, (ssp, elem)]
            cum_TEs.loc[year:, (ssp, perc)] = cum_TEs.loc[year:, (ssp, perc)] + 1

# %% calculate slow tipping probabilities 
P_sample = pickle.load(open('Data/Params/coupled_ensemble/P_sample.pkl', 'rb'))

for ssp in tip_years.columns.levels[0].values:
    for elem in elements:
        T_100 = T[ssp].loc[2400:].mean(axis=0) 
        p_tip = ((T_100 > P_sample[elem]) * 1).mean()
        slow_tip_probs.loc[:, (ssp, elem)] = p_tip

# %% save tipping probs and years
tip_years.to_csv('Data/model_output/' + ens_name + '/Y_tip.csv')
tip_emms.to_csv('Data/model_output/' + ens_name + '/Cum_C_tip.csv')
cum_TEs.to_csv('Data/model_output/' + ens_name + '/Cum_TEs.csv')
slow_tip_probs.to_csv('Data/model_output/' + ens_name + '/slow_tip_probs.csv')


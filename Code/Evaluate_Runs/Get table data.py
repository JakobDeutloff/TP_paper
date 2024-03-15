"""
Collect data for Table 2
"""
# %%
from Code.Read_data.Read_SSP_output import read_full_T, read_SSP_outout, read_TE_output
import pandas as pd


# %% load data
res_tip = read_SSP_outout('coupled_ensemble')
T_full_tip = read_full_T('coupled_ensemble')
T_full_const = read_full_T('uncoupled_ensemble')
TEs = read_TE_output('coupled_ensemble')

# %% Get T data for table

# Calculate differences of T
diff_T = T_full_tip - T_full_const

# Set up table DF
years = [2050, 2100, 2200, 2300, 2400, 2500]
ssps = ['ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585']
vals = ['T', 'carbon_dioxide', 'methane']
cols = pd.MultiIndex.from_product([ssps, vals, years])
quants = [0.05, 0.5, 0.95]
table_data = pd.DataFrame(index=quants, columns=cols)

# Set up df for cumulative emms
cum_CO2 = {}
cum_CH4 = {}

# %% calculate cumulative emissions
for ssp in ssps:
    cum_CO2[ssp] = \
        TEs['TE_total'].xs((ssp, 'CO2_emm'), axis=1, level=(0, 3)).groupby(level=0, axis=1).sum().cumsum()
    cum_CH4[ssp] = \
        TEs['TE_total'].xs((ssp, 'CH4_emm'), axis=1, level=(0, 3)).groupby(level=0, axis=1).sum().cumsum()

# %% fill table
for ssp in ssps:
    for year in years:
        for quant in quants:
            table_data.loc[quant, (ssp, 'T', year)] = diff_T.loc[year, ssp].quantile(quant).round(2)
            table_data.loc[quant, (ssp, 'carbon_dioxide', year)] = cum_CO2[ssp].loc[year].quantile(quant).round(0)
            table_data.loc[quant, (ssp, 'methane', year)] = (cum_CH4[ssp].loc[year].quantile(quant) * 1e-3).round(1)



# %%

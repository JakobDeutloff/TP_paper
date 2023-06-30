from Code.Read_data.Read_SSP_output import read_full_T, read_SSP_outout, read_TE_output
import pandas as pd
# %% load data
res_tip = read_SSP_outout('5000_tip_2')
T_full_tip = read_full_T('5000_tip_2')
T_full_const = read_full_T('5000_const_2')
TEs = read_TE_output('5000_tip_2')

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

# %% check CH4 fractions

def CH4_to_C(CH4):
    """
    :param CH4: in MtCH4
    :return: C in GtC
    """
    weight_C = 12.011  # weight of carbon in u
    weight_CH4 = 16  # weight of CH4 in u
    factor = (weight_CH4 / weight_C) * 1e3
    C = CH4 / factor
    return C

def get_frac(CH4, CO2):

    CH4_C = CH4_to_C(CH4*1e3)
    print(CH4_C/(CO2 + CH4_C))



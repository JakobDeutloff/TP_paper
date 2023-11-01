"""
Calculate carbon emissions from CO2 and CH4 emissions, plot carbon emissions (Fig. 3, Fig. S9),
plot atmospheric concentrations of CO2 and CH4 (Fig. S10), calculate emission shares used in the text.
"""
# %%
import matplotlib.pyplot as plt
import matplotlib as mpl
from Code.Read_data.Read_SSP_output import read_SSP_outout, read_TE_output
import pandas as pd
import numpy as np
# mpl.use('Qt5Agg')  # Use interactive plotting interface


# %%
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


# %% Read data
const = read_SSP_outout('uncoupled_ensemble')
tip = read_SSP_outout('coupled_ensemble')
TEs = read_TE_output('coupled_ensemble')

# %% Calculate carbon emissions

# Only look at tier 1 SSPs
ssps = ['ssp126', 'ssp245', 'ssp370', 'ssp585']
elems = ['AMAZ', 'PFTP', 'PFAT']
percentiles = ['0.05', '0.5', '0.95']

# Initialize dicts to store data for table
multiidx = pd.MultiIndex.from_product([ssps, ['CO2', 'CH4']])
emm_const = pd.DataFrame(index=const['Emm'].index, columns=multiidx)
emm_const.loc[:, :] = 0

multiidx = pd.MultiIndex.from_product([ssps, percentiles, ['CO2', 'CH4']])
emm_tip = pd.DataFrame(index=const['Emm'].index, columns=multiidx)
emm_tip.loc[:, :] = 0

multiidx = pd.MultiIndex.from_product([ssps, percentiles])
cum_emm_tip = pd.DataFrame(index=const['Emm'].index, columns=multiidx)
cum_emm_tip.loc[:, :] = 0

# Calculate emissions

for i in range(4):

    # Uncoupled
    emm_const.loc[:, (ssps[i], 'CO2')] = const['Emm'].loc[:, (ssps[i], 'carbon_dioxide')]
    emm_const.loc[:, (ssps[i], 'CH4')] = CH4_to_C(const['Emm'].loc[:, (ssps[i], 'methane')])

    # Coupled
    for perc in percentiles:
        for elem in elems:
            # Emm
            emm_tip.loc[:, (ssps[i], perc, 'CO2')] += TEs['TE'].loc[:, (ssps[i], elem, 'CO2_emm', perc)]
            emm_tip.loc[:, (ssps[i], perc, 'CH4')] += CH4_to_C(TEs['TE'].loc[:, (ssps[i], elem, 'CH4_emm', perc)])

            # Cumulative Emm
            cum_emm_tip.loc[:, (ssps[i], perc)] += TEs['TE'].loc[:, (ssps[i], elem, 'C_stock', perc)]

# %% Plot C emissions together with SSP emissions
fig, axes = plt.subplots(1, 4, figsize=(13, 4), sharex='col', sharey='row')
start = 2000
end = 2500
ssps_plot = ['ssp126', 'ssp245', 'ssp370', 'ssp585']

for i in range(4):
    # Plot Cumulative C Emissions
    ax = axes[i]
    # Uncoupled
    ax.plot(emm_const.loc[start:end, ssps_plot[i]].sum(axis=1).cumsum(), color='k', linestyle='--',
            label='SSP Emissions')
    # Coupled
    ax.plot(cum_emm_tip.loc[start:end, (ssps_plot[i], '0.5')], color='red', label='Median CTEM Emissions')
    ax.fill_between(cum_emm_tip.loc[start:end].index,
                    cum_emm_tip.loc[start:end, (ssps_plot[i], '0.05')].astype('float'),
                    cum_emm_tip.loc[start:end, (ssps_plot[i], '0.95')].astype('float'),
                    color='red', alpha=0.5, linestyle='-', label='5th to 95th percentile of CTEM Emissions')
    ax.set_xlabel('Year')
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.set_title(ssps_plot[i])

axes[0].set_ylabel('Cumulative C Emissions [GtC]')

# Legend
handles, labels = axes[1].get_legend_handles_labels()
# fig.subplots_adjust(bottom=0.2)
plt.legend(handles, labels, ncol=3, bbox_to_anchor=(0.2, -0.15))
plt.savefig('Plots/Carbon_emms/C_emm_ssps.png', dpi=300, bbox_inches='tight')
plt.show()


# %% Plot CO2 and CH4 concentrations
fig, axes = plt.subplots(2, 4, figsize=(13, 7), sharex='col')
start = 2000
end = 2500
ssps_plot = ['ssp126', 'ssp245', 'ssp370', 'ssp585']
for i in range(4):
    # Plot CO2 anomaly
    ax = axes[0, i]
    CO2_05_c = const['C'].loc[start:end, (ssps_plot[i], 'carbon_dioxide', '0.05')]
    CO2_50_c = const['C'].loc[start:end, (ssps_plot[i], 'carbon_dioxide', '0.5')]
    CO2_95_c = const['C'].loc[start:end, (ssps_plot[i], 'carbon_dioxide', '0.95')]
    CO2_05 = tip['C'].loc[start:end, (ssps_plot[i], 'carbon_dioxide', '0.05')]
    CO2_50 = tip['C'].loc[start:end, (ssps_plot[i], 'carbon_dioxide', '0.5')]
    CO2_95 = tip['C'].loc[start:end, (ssps_plot[i], 'carbon_dioxide', '0.95')]

    ax.plot(CO2_50, color='r', label='Median Coupled')
    ax.fill_between(CO2_50.index, CO2_05, CO2_95, color='r', alpha=0.3, linestyle='-',
                    label='5th to 95th Percentile Coupled')
    ax.plot(CO2_50_c, color='k', linestyle='--', label='Median Uncoupled')
    ax.fill_between(CO2_50.index, CO2_05_c, CO2_95_c, color='grey', alpha=0.3, linestyle='-',
                    label='5th to 95th Percentile Uncoupled')
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.set_title(ssps_plot[i])

    # Plot CH4 anomaly
    ax = axes[1, i]
    CH4_05_c = const['C'].loc[start:end, (ssps_plot[i], 'methane', '0.05')]
    CH4_50_c = const['C'].loc[start:end, (ssps_plot[i], 'methane', '0.5')]
    CH4_95_c = const['C'].loc[start:end, (ssps_plot[i], 'methane', '0.95')]
    CH4_05 = tip['C'].loc[start:end, (ssps_plot[i], 'methane', '0.05')]
    CH4_50 = tip['C'].loc[start:end, (ssps_plot[i], 'methane', '0.5')]
    CH4_95 = tip['C'].loc[start:end, (ssps_plot[i], 'methane', '0.95')]

    ax.plot(CH4_50, color='r')
    ax.fill_between(CH4_50.index, CH4_05, CH4_95, color='r', alpha=0.3, linestyle='-')
    ax.plot(CH4_50_c, color='k', linestyle='--')
    ax.fill_between(CH4_50.index, CH4_05_c, CH4_95_c, color='grey', alpha=0.3, linestyle='-')
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.set_xlabel('Year')

handles, labels = axes[0, 0].get_legend_handles_labels()
plt.legend(handles, labels, ncol=4, bbox_to_anchor=(0.8, -0.19))

axes[0, 0].set_ylabel(r'$\mathrm{CO_2}$ Concentration [ppm]')
axes[1, 0].set_ylabel(r'$\mathrm{CH_4}$ Concentration [ppb]')

fig.savefig('Plots/Carbon_emms/concentrations.png', dpi=300, bbox_inches='tight')
fig.show()


# %% Share of Cumulative Emissions
start = 2000
end = 2500

max = cum_emm_tip.loc[:, ('ssp126', '0.95')].max()
const = emm_const.loc[start:end, ('ssp126')].sum(axis=1).cumsum().iloc[-1]
share = max / const
print('share: ' + str(share) + ', max: ' + str(max))

max = cum_emm_tip.loc[:, ('ssp245', '0.95')].max()
const = emm_const.loc[start:end, ('ssp245')].sum(axis=1).cumsum().iloc[-1]
share = max / const
print('share: ' + str(share) + ', max: ' + str(max))

max = cum_emm_tip.loc[:, ('ssp370', '0.95')].max()
const = emm_const.loc[start:end, ('ssp370')].sum(axis=1).cumsum().iloc[-1]
share = max / const
print('share: ' + str(share) + ', max: ' + str(max))

max = cum_emm_tip.loc[:, ('ssp585', '0.95')].max()
const = emm_const.loc[start:end, ('ssp585')].sum(axis=1).cumsum().iloc[-1]
share = max / const
print('share: ' + str(share) + ', max: ' + str(max))


# %%

import matplotlib.pyplot as plt
from Code.Read_data.Read_SSP_output import read_SSP_outout, read_rand_mems, read_TE_output
import numpy as np
import pandas as pd
import matplotlib as mpl

# %% read data
const = read_SSP_outout('5000_const_2')
tip = read_SSP_outout('5000_tip_2')
rand_mems = read_rand_mems('5000_tip_2')
rand_mems_const = read_rand_mems('5000_const_2')
TE_out = read_TE_output('5000_tip_2')


# %% Define Plot functions
# Function for T comparison
def plot_T(ssp, ax):
    # historical
    ax.plot(const['T'].loc[1850:2016, ('ssp245', '0.5')], color='grey', label='Median Historical')
    ax.fill_between(np.arange(1850, 2017), const['T'].loc[1850:2016, ('ssp245', '0.05')].values,
                    const['T'].loc[1850:2016, ('ssp245', '0.95')].values, color='grey', alpha=0.3, linestyle='-',
                    label='5th to 95th Percentile Historical')
    # Uncoupled
    ax.plot(const['T'].loc[2016:, (ssp, '0.5')], label='Median Uncoupled', color='k', linestyle='--')
    ax.fill_between(const['T'].loc[2016:].index, const['T'].loc[2016:, (ssp, '0.05')].values,
                    const['T'].loc[2016:, (ssp, '0.95')].values, color='k', alpha=0.3, linestyle='-',
                    label='5th to 95th Percentile Uncoupled')
    # Coupled
    ax.plot(tip['T'].loc[2016:, (ssp, '0.5')], label='Median Coupled', color='r')
    ax.fill_between(np.arange(2016, 2501), tip['T'].loc[2016:, (ssp, '0.05')].values,
                    tip['T'].loc[2016:, (ssp, '0.95')].values, color='r', alpha=0.3, linestyle='-',
                    label='5th to 95th Percentile Coupled')

    ax.set_title(ssp)


def plot_TvsE(ssp, ax):
    end_year = 2250
    colors = ['red', 'orange', 'green', 'yellow', 'k', 'purple', 'deeppink', 'magenta', 'cyan', 'blue']
    mems = np.array([2, 4, 5, 11, 30])

    # Calculate cum Emissions
    weight_C = 12.011  # weight of carbon in u
    weight_CH4 = 16  # weight of CH4 in u
    fac_CH4 = (weight_C / weight_CH4) * 1e3  # used in fair for conversion from GtC to MtCH4
    em_const = const['Emm'][ssp].loc[:end_year, 'carbon_dioxide'] # + const['Emm'][ssp].loc[:end_year, 'methane'] * (1 / fac_CH4)
    cum_em_const = em_const.cumsum()
    cum_emms_SSP[ssp] = cum_em_const

    # Hist
    ax.plot(cum_em_const.loc[:2016], const['T'][ssp].loc[:2016, '0.5'], color='grey', label='Median Historical')
    ax.fill_between(cum_em_const.loc[:2016], const['T'][ssp].loc[:2016, '0.05'],
                    const['T'][ssp].loc[:2016, '0.95'],
                    color='grey', alpha=0.3, linestyle='-', label='5th to 95th Percentile Historical')

    # Const
    ax.plot(cum_em_const.loc[2016:], const['T'][ssp].loc[2016:end_year, '0.5'], color='k', linestyle='--',
            label='Median Uncoupled')
    ax.fill_between(cum_em_const.loc[2016:], const['T'][ssp].loc[2016:end_year, '0.05'],
                    const['T'][ssp].loc[2016:end_year, '0.95'],
                    color='k', alpha=0.3, linestyle='-', label='5th to 95th Percentile Uncoupled')
    # Tip
    ax.plot(cum_em_const.loc[2016:], tip['T'][ssp].loc[2016:end_year, '0.5'], color='r', label='Median Coupled')
    ax.fill_between(cum_em_const.loc[2016:], tip['T'][ssp].loc[2016:end_year, '0.05'],
                    tip['T'][ssp].loc[2016:end_year, '0.95'],
                    color='r', alpha=0.3, linestyle='-', label='5th to 95th Percentile Coupled')

    # Rand Members
    for i in range(len(mems)):
        ax.plot(cum_em_const.loc[2016:], rand_mems[ssp]['T'].loc[2016:end_year].iloc[:, mems[i]], color='red',
                linewidth=0.5, linestyle='--')
        ax.plot(cum_em_const.loc[2016:], rand_mems_const[ssp]['T'].loc[2016:end_year].iloc[:, mems[i]], color='k',
                linewidth=0.5, linestyle='--')

    ax.set_xlabel(r'Cumulative $\mathrm{CO_2}$ Emissions [GtC]')
    ax.set_title(ssp)

    return cum_emms_SSP


# %% Plot T comparison
mpl.use('Qt5Agg')
ssps = ['ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585']
fig, axes = plt.subplots(2, 3, figsize=(11, 7), sharex='col')
for i in range(5):
    plot_T(ssps[i], axes.flatten()[i])

axes[1, 2].remove()
axes[0, 2].tick_params(labelbottom=True)
axes[0, 0].set_ylabel('GMST anomaly [°C]')
axes[1, 0].set_ylabel('GMST anomaly [°C]')
axes[1, 0].set_xlabel('Year')
axes[1, 1].set_xlabel('Year')
axes[0, 2].set_xlabel('Year')

handles, labels = axes[1, 1].get_legend_handles_labels()
fig.legend(handles, labels, ncol=1, bbox_to_anchor=(0.97, 0.47))
fig.tight_layout()
fig.savefig('Plots/Comparison/5000_tip_2_5000_const_2/Temp/Tier1_temp.png', dpi=300)
plt.show()

# %% Plot Tvs E
ssps = ['ssp245', 'ssp370', 'ssp585']
mpl.use('Qt5Agg')
# fake ens mems entries for legend
fig_2, ax = plt.subplots()
ax.plot([1, 2, 3], [1, 2, 3], linewidth=0.4, linestyle='--', color='k', label='Uncoupled Ens. Mem.')
ax.plot([1, 2, 3], [1, 2, 3], linewidth=0.4, linestyle='--', color='r', label='Coupled Ens. Mem.')
handles_2, labels_2 = ax.get_legend_handles_labels()
fig_2.clf()
ax.cla()

# Plot T vs E
fig, axes = plt.subplots(1, 3, figsize=(11, 4))
cum_emms_SSP = {}
for i in range(3):
    plot_TvsE(ssps[i], axes[i])

axes[0].set_ylabel('GMST anomaly [°C]')
# Legend
handles, labels = axes[1].get_legend_handles_labels()
handles = handles + handles_2
labels = labels + labels_2
fig.subplots_adjust(bottom=0.2)
plt.legend(handles, labels, ncol=3, bbox_to_anchor=(0.8, -0.17))

fig.savefig('Plots/Comparison/5000_tip_2_5000_const_2/Temp/TvsE.png', dpi=300, bbox_inches='tight')
fig.show()


# %% Collect data for T table


# Set up dataframe
SSPs = ['ssp119', 'ssp245', 'ssp585']
ens = ['coupled', 'uncoupled']
quants = ['0.5', '0.95', 'PFTP', 'PFAT', 'AMAZ']
years = [2050, 2100, 2200, 2300, 2400, 2500]
columns = pd.MultiIndex.from_product([SSPs, ens, quants])
T_vals = pd.DataFrame(index=years, columns=columns)

# Collect values
for ssp in SSPs:
    for year in years:
        for perc in ['0.5', '0.95']:
            T_vals.loc[year, (ssp, 'coupled', perc)] = tip['T'].loc[year, (ssp, perc)]
            T_vals.loc[year, (ssp, 'uncoupled', perc)] = const['T'].loc[year, (ssp, perc)]
        for elem in ['PFTP', 'PFAT', 'AMAZ']:
            T_vals.loc[year, (ssp, 'coupled', elem)] = \
                (TE_out['TE_total'].xs((ssp, elem, 'bool_tip'), axis=1, level=(0, 2, 3)).loc[year].sum() / 5000) * 100

# calculate differences
diff_5 = T_vals.xs(('coupled', '0.5'), axis=1, level=(1, 2)) - T_vals.xs(('uncoupled', '0.5'), axis=1, level=(1, 2))
diff_95 = T_vals.xs(('coupled', '0.95'), axis=1, level=(1, 2)) - T_vals.xs(('uncoupled', '0.95'), axis=1, level=(1, 2))
perc_5 = (diff_5 / T_vals.xs(('uncoupled', '0.5'), axis=1, level=(1, 2))) * 100
perc_95 = (diff_95 / T_vals.xs(('uncoupled', '0.95'), axis=1, level=(1, 2))) * 100

# %% Investigate which random ensemble members in TvsE feature tipping

mems = np.array([2, 4, 5, 11, 30])  #np.array([2, 12, 17, 11, 30])
ssp = 'ssp245'
res = {'AMAZ': np.zeros(5), 'PFTP': np.zeros(5), 'PFAT': np.zeros(5)}

members = list(rand_mems[ssp]['tipping_elements'].xs(('PFTP', 'bool_tip'), axis=1, level=(1, 2)).columns)
for elem in ['PFTP', 'AMAZ', 'PFAT']:
    for i in range(len(mems)):
        res[elem][i] = rand_mems[ssp]['tipping_elements'].xs((members[mems[i]], elem, 'bool_tip'), axis=1,
                                                             level=(0, 1, 2)).max()

    print(elem + ': ' + str(res[elem]))

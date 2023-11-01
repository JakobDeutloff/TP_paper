"""
Plot temperature from the coupled and the uncoupled ensemble (Fig. 4)
"""
# %%
import matplotlib.pyplot as plt
from Code.Read_data.Read_SSP_output import read_SSP_outout, read_rand_mems, read_TE_output
import numpy as np
import pandas as pd
import matplotlib as mpl

# %% read data
const = read_SSP_outout('uncoupled_ensemble')
tip = read_SSP_outout('coupled_ensemble')
rand_mems = read_rand_mems('coupled_ensemble')
rand_mems_const = read_rand_mems('uncoupled_ensemble')
TE_out = read_TE_output('coupled_ensemble')


#%% Function for T comparison
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


# %% Plot T comparison
# mpl.use('Qt5Agg')
ssps = ['ssp126', 'ssp245', 'ssp370', 'ssp585']
fig, axes = plt.subplots(1, 4, figsize=(13, 4), sharex='col', sharey='row')
for i in range(4):
    plot_T(ssps[i], axes[i])

axes[0].set_ylabel('GMST anomaly [°C]')
for ax in axes:
    ax.set_xlabel('Year')

# make legend at the bottom of the figure
handles, labels = axes[1].get_legend_handles_labels()
fig.legend(handles, labels, ncol=3, bbox_to_anchor=(0.825, 0))
fig.savefig('Plots/Temp/Tier1_temp.png', dpi=300, bbox_inches='tight')
plt.show()





# %%

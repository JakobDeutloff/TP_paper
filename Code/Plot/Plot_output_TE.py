"""
Plot Output crom CTEM - C Stock, CO2 emms and CH4 emms
"""

from Code.Read_data.Read_SSP_output import read_TE_output, read_SSP_outout
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

# %% Load data
ens_name = '5000_tip_2'
N_sets = 5000
TE_out = read_TE_output(ens_name)
SSP_out = read_SSP_outout('5000_const')
TE_total = TE_out['TE_total']


# %% define functions

def plot_ens_mem(TE_out, ax):
    start = 2000
    TE_out.quantile(0.5, axis=1).loc[start:].plot(ax=ax, color='k', label='Median')
    ax.fill_between(TE_out.loc[start:].index, TE_out.quantile(0.05, axis=1).loc[start:],
                    TE_out.quantile(0.95, axis=1).loc[start:], color='k', alpha=0.4,
                    label='5th to 95th Percentile')
    ax.plot(TE_out.iloc[:, random_members].loc[start:], color='r', linewidth=0.5, label='Ensemble Member')


def plot_ssp(ssp):
    fig, axes = plt.subplots(3, 3, figsize=(10, 7), sharex='col')

    # C Stock
    plot_ens_mem(TE_total.xs((ssp, 'PFTP', 'C_stock'), level=(0, 2, 3), axis=1), axes[0, 0])
    plot_ens_mem(TE_total.xs((ssp, 'AMAZ', 'C_stock'), level=(0, 2, 3), axis=1), axes[0, 1])
    plot_ens_mem(TE_total.xs((ssp, 'PFAT', 'C_stock'), level=(0, 2, 3), axis=1), axes[0, 2])
    axes[0, 0].set_ylabel('Cumulative C Em. [GtC]')

    # CO2 emm
    plot_ens_mem(TE_total.xs((ssp, 'PFTP', 'CO2_emm'), level=(0, 2, 3), axis=1), axes[1, 0])
    plot_ens_mem(TE_total.xs((ssp, 'AMAZ', 'CO2_emm'), level=(0, 2, 3), axis=1), axes[1, 1])
    plot_ens_mem(TE_total.xs((ssp, 'PFAT', 'CO2_emm'), level=(0, 2, 3), axis=1), axes[1, 2])
    axes[1, 0].set_ylabel(r'$\mathrm{CO_2}$ Em. [$\mathrm{GtC ~ yr^{-1}}$]')

    # CH4 emm
    plot_ens_mem(TE_total.xs((ssp, 'PFTP', 'CH4_emm'), level=(0, 2, 3), axis=1), axes[2, 0])
    plot_ens_mem(TE_total.xs((ssp, 'AMAZ', 'CH4_emm'), level=(0, 2, 3), axis=1), axes[2, 1])
    plot_ens_mem(TE_total.xs((ssp, 'PFAT', 'CH4_emm'), level=(0, 2, 3), axis=1), axes[2, 2])
    axes[2, 0].set_ylabel(r'$\mathrm{CH_4}$ Em. [$\mathrm{MtCH_4 ~ yr^{-1}}$]')

    # Labeling
    axes[0, 0].set_title('PFTP')
    axes[0, 1].set_title('AMAZ')
    axes[0, 2].set_title('PFAT')
    axes[2, 0].set_xlabel('Year')
    axes[2, 1].remove()
    axes[2, 2].set_xlabel('Year')
    axes[1, 1].tick_params(labelbottom=True)
    axes[1, 1].set_xlabel('Year', visible='true')

    handles, labels = axes[1, 1].get_legend_handles_labels()
    fig.legend(handles[:3], labels[:3], ncol=1, bbox_to_anchor=(0.625, 0.3))
    fig.tight_layout()
    plt.savefig('Plots/SSP_ensembles/' + ens_name + '/CTEM/CTEM_' + ssp + '.png', dpi=300)
    plt.show()
    # fig.clf()


def plot_ssp119():
    fig, axes = plt.subplots(3, 1, figsize=(5.5, 7), sharex='col')
    ssp = 'ssp119'

    # C Stock
    plot_ens_mem(TE_total.xs((ssp, 'PFAT', 'C_stock'), level=(0, 2, 3), axis=1), axes[0])
    axes[0].set_ylabel('Cumulative C Em. [GtC]')

    # CO2 emm
    plot_ens_mem(TE_total.xs((ssp, 'PFAT', 'CO2_emm'), level=(0, 2, 3), axis=1), axes[1])
    axes[1].set_ylabel(r'$\mathrm{CO_2}$ Em. [$\mathrm{GtC ~ yr^{-1}}$]')

    # CH4 emm
    plot_ens_mem(TE_total.xs((ssp, 'PFAT', 'CH4_emm'), level=(0, 2, 3), axis=1), axes[2])
    axes[2].set_ylabel(r'$\mathrm{CH_4}$ Em. [$ \mathrm{MtCH_4 ~ yr^{-1}}$]')

    # Labeling
    axes[0].set_title('PFAT')

    handles, labels = axes[0].get_legend_handles_labels()
    plt.legend(handles[:3], labels[:3], ncol=1, bbox_to_anchor=(0.8, -0.25))
    plt.savefig('Plots/SSP_ensembles/' + ens_name + '/CTEM/CTEM_' + ssp + '2.png', dpi=300, bbox_inches='tight')
    fig.show()
    # fig.clf()


# %% Plot all SSPs

# Choose 10 random ensemble members - same for every ssp
random_members = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
mpl.use('Qt5Agg')  # Use interactive plotting interface
for ssp in ['ssp585', 'ssp245']: #list(TE_total.columns.levels[0]):
    plot_ssp(ssp)
plot_ssp119()


# %% Analyze values
def give_tipped(ssp):
    for elem in ['PFAT', 'AMAZ', 'PFTP']:
        num = (TE_total.xs((ssp, elem, 'bool_tip'), level=(0, 2, 3), axis=1).sum(axis=1).iloc[-1] / 5000) * 100
        print(ssp + ' ' + elem + ' tipped : ' + str(num))


for ssp in ['ssp119', 'ssp245', 'ssp585']:
    give_tipped(ssp)


# %% Calculate maximum emission peaks under ssp585

def max_emms_mag(elem):
    # Get emission peaks
    CO2 = TE_total.xs(('ssp585', elem, 'CO2_emm'), level=(0, 2, 3), axis=1).max()
    CO2 = CO2[CO2 > 0]
    CH4 = TE_total.xs(('ssp585', elem, 'CH4_emm'), level=(0, 2, 3), axis=1).max()
    CH4 = CH4[CH4 > 0]

    # Get index of emm peaks
    CO2_idx = []
    for mem in CO2.index.values:
        vals = TE_total.xs(('ssp585', mem, elem, 'CO2_emm'), level=(0, 1, 2, 3), axis=1)
        CO2_idx += [vals[(CO2[mem] == vals).values].index.values[0]]

    # Exclude members with max emms in 2500
    CO2 = CO2[np.array(CO2_idx) != 2500]
    CO2_idx = np.array(CO2_idx)[np.array(CO2_idx) != 2500]

    # Magnitude of emm peaks
    print('SSP585 magnitudes of ' + elem + ': CO2: ' + str(np.mean(CO2)) + ' (' + str(np.min(CO2)) + ', '
          + str(np.max(CO2)) + ')' + ' ,CH4: ' + str(np.mean(CH4)) + ' (' + str(np.min(CH4)) + ', '
          + str(np.max(CH4)) + ')')
    # Timing of emm peaks
    print('SSP585 timing of ' + elem + ': ' + str(np.mean(CO2_idx)) + ' (' + str(np.min(CO2_idx)) + ', '
          + str(np.max(CO2_idx)) + ')')


max_emms_mag('PFTP')
max_emms_mag('AMAZ')
max_emms_mag('PFAT')

# %% Get anthropogenic emission peaks

CO2_ssp = SSP_out['Emm'].loc[:, ('ssp585', 'carbon_dioxide')]
CH4_ssp = SSP_out['Emm'].loc[:, ('ssp585', 'methane')]
CO2_max = CO2_ssp.max()
CH4_max = CH4_ssp.max()
CO2_idx = CO2_ssp[CO2_ssp == CO2_max].index.values[0]
CH4_idx = CH4_ssp[CH4_ssp == CH4_max].index.values[0]

print('SSP emissions: CO2: ' + str(CO2_max) + ' in ' + str(CO2_idx) + ' , CH4: ' + str(CH4_max) + ' in ' + str(CH4_idx))

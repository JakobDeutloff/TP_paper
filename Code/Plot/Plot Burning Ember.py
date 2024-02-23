"""
Create Burning Ember plots (Fig. 6, 7, 8, S11 - S17) and calculate how much earlier tipping happens on average due to carbon TEs
"""

# %%
from Code.Read_data.Read_SSP_output import read_SSP_outout, read_probs
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import pandas as pd

# %% Read data
ens_name_tip = 'coupled_ensemble'
ens_name_const = 'uncoupled_ensemble'

SSP_out_tip = read_SSP_outout(ens_name_tip)
Probs_tip = read_probs(ens_name_tip)
SSP_out_const = read_SSP_outout(ens_name_const)
Probs_const = read_probs(ens_name_const)


# %% Define Functions

def plot_burning_emeber(ssp, max_year=2500):
    columns = Probs_const['Years'].columns.levels[1]
    values_max_const = SSP_out_const['Tip_probs_total'].loc[:max_year][ssp][columns].max().values
    values_max_tip = SSP_out_tip['Tip_probs_total'].loc[:max_year][ssp][columns].max().values
    values_slow_tipping = Probs_const['slow_tip_probs'].loc[2450, ssp][columns].values

    fig, ax = plt.subplots(figsize=(8, 5))

    # set grid
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.05, 0.33, 0.5, 0.66, 0.95, 1])
    ax.grid(which='major', axis='y', linestyle='--', color='gray', zorder=0)

    # Plot bars of coupled run
    ax.bar(columns, values_max_tip, color='red', edgecolor='grey', zorder=2)
    # plot instantaneous tipping probabilities (uncoupled)
    ax.bar(columns, values_max_const, color="orange", edgecolor='grey', zorder=3)
    # plot slow tipping probabilities (uncoupled)
    if max_year == 2500:
        ax.bar(columns, values_slow_tipping, color="gold", edgecolor='grey', zorder=4)
    first, last = ax.get_xlim()
    plt.xlim(left=first + 0.7, right=last - 0.7)
    ax.set_ylabel('Probability of Tipping')

    ax.set_title(ssp)
    plt.tight_layout()

# %% Plot Burning Embers diagramm

for ssp in Probs_const['Years'].columns.levels[0]:
    plot_burning_emeber(ssp, 2200)
    plt.tight_layout()
    plt.savefig('Plots/Burning_ember_2200/' + ssp + '_200.png', dpi=300)

# %% calculate total increase in probability of tippig for each SSP
elems = list(Probs_const['Years'].columns.levels[1])
ssps = list(Probs_const['Years'].columns.levels[0])
multicols = pd.MultiIndex.from_product([ssps, elems])

P_const = pd.DataFrame(index=[0], columns=multicols)
P_const_200 = pd.DataFrame(index=[0], columns=multicols)
P_tip = pd.DataFrame(index=[0], columns=multicols)
P_mean_diff = pd.DataFrame(index=[0], columns=ssps)
P_mean_const = pd.DataFrame(index=[0], columns=ssps)
P_diff = pd.DataFrame(index=[0], columns=multicols)

for ssp in ssps:
    for elem in elems:
        P_const.loc[:, (ssp, elem)] = SSP_out_const['Tip_probs_total'][ssp][elem].max()
        P_const_200.loc[:, (ssp, elem)] = SSP_out_const['Tip_probs_total'][ssp][elem].loc[:2200].max()
        P_tip.loc[:, (ssp, elem)] = SSP_out_tip['Tip_probs_total'][ssp][elem].max()
        P_diff.loc[:, (ssp, elem)] = P_tip.loc[:, (ssp, elem)] - P_const.loc[:, (ssp, elem)]

    P_mean_diff.loc[:, ssp] = P_tip[ssp].mean(axis=1) - P_const[ssp].mean(axis=1)
    P_mean_const.loc[:, ssp] = P_const[ssp].mean(axis=1)


# %% calculate how much earlier tipping happens on average

Y_diff_5 = {}
Y_diff_95 = {}


for ssp in Probs_const['Years'].columns.levels[0]:
    yrs_5 = []
    yrs_95 = []
    for elem in Probs_const['Years'].columns.levels[1]:
        yrs_5 += [Probs_tip['Years'].loc[0.5, (ssp, elem)] - Probs_const['Years'].loc[0.5, (ssp, elem)]]
        yrs_95 += [Probs_tip['Years'].loc[0.95, (ssp, elem)] - Probs_const['Years'].loc[0.95, (ssp, elem)]]
        if yrs_5[-1] < -10:
            print(ssp + ' ' + elem + ' ' + str(yrs_5[-1]))

    Y_diff_5[ssp] = yrs_5
    Y_diff_95[ssp] = yrs_95
# %%
print('Mean Diff [years] at 50%: ' + str(np.nanmean(list(Y_diff_5.values()))))
print('Mean Diff [years] at 95%: ' + str(np.nanmean(list(Y_diff_95.values()))))

# %%

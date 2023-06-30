"""
Make comparison plots fron tipping points and normal run. Burning ember and temperature
"""

from Code.Read_data.Read_SSP_output import read_SSP_outout, read_probs
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import pandas as pd

# %% Read data
ens_name_tip = '5000_tip_2'
ens_name_const = '5000_const_2'

SSP_out_tip = read_SSP_outout(ens_name_tip)
Probs_tip = read_probs(ens_name_tip)
SSP_out_const = read_SSP_outout(ens_name_const)
Probs_const = read_probs(ens_name_const)


# %% Define Functions

def gradientbars(bars, cmap):
    ax = bars[0].axes
    lim = ax.get_xlim() + ax.get_ylim()
    ax.axis(lim)
    for bar in bars:
        bar.set_facecolor("none")
        x, y = bar.get_xy()
        w, h = bar.get_width(), bar.get_height()
        grad = np.atleast_2d(np.linspace(0, 1 * h / 1, 1000)).T
        # zorder of 2 to get gradients above the facecolor, but below the bar outlines
        ax.imshow(grad, extent=[x, x + w, y, y + h], origin='lower', aspect="auto", zorder=2,
                  norm=cm.colors.NoNorm(vmin=0, vmax=1), cmap=plt.get_cmap(cmap))


def plot_burning_emeber(ssp, max_year=2500):
    columns = Probs_const['Years'].columns.levels[1]
    rows = ['0.5', '0.5', '0.95', '0.95']
    values_max_const = SSP_out_const['Tip_probs_total'].loc[:max_year][ssp][columns].max().values
    values_max_tip = SSP_out_tip['Tip_probs_total'].loc[:max_year][ssp][columns].max().values

    fig, ax = plt.subplots(figsize=(8, 5))

    # set grid
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.05, 0.33, 0.5, 0.66, 0.95, 1])
    ax.grid(which='major', axis='y', linestyle='--', color='gray', zorder=0)

    # Plot bars of coupled run
    ax.bar(columns, values_max_tip, color='mediumvioletred', edgecolor='grey', zorder=2)
    # plot bars of standard run
    my_bar = ax.bar(columns, values_max_const, color='grey', edgecolor='grey', zorder=3)

    gradientbars(my_bar, 'YlOrRd')
    plt.xticks([])
    first, last = ax.get_xlim()
    plt.xlim(left=first + 0.7, right=last - 0.7)
    ax.set_ylabel('Probability of Tipping')

    # draw table
    # modify table text
    text_const = Probs_const['Years'][ssp][columns].loc[(0.5, 0.95), :].values
    text_tip = Probs_tip['Years'][ssp][columns].loc[(0.5, 0.95), :].values
    # Exclude values beyond max_year
    text_const[text_const > max_year] = np.nan
    text_tip[text_tip > max_year] = np.nan

    # construct array of strings or ompty strings for nan
    text_str_const = [[None for _ in range(text_const.shape[1])] for _ in range(text_const.shape[0])]
    for i in range(text_const.shape[0]):
        for j in range(text_const.shape[1]):
            try:
                text_str_const[i][j] = str(int(text_const[i, j]))
            except:
                text_str_const[i][j] = ''

    text_str_tip = [[None for _ in range(text_tip.shape[1])] for _ in range(text_tip.shape[0])]
    for i in range(text_tip.shape[0]):
        for j in range(text_tip.shape[1]):
            try:
                text_str_tip[i][j] = str(int(text_tip[i, j]))
            except:
                text_str_tip[i][j] = ''

    # Merge strs
    str_table = [[None for _ in range(text_tip.shape[1])] for _ in range(text_tip.shape[0]*2)]
    for i in [0, 2]:
        str_table[i] = text_str_tip[int(i/2)]
        str_table[i+1] = text_str_const[int(i/2)]

    the_table = ax.table(cellText=str_table, rowLabels=rows, colLabels=columns,
                         loc='bottom', cellLoc='center', colLoc='center', rowLoc='center', edges='horizontal')

    # Set cell heights
    cellDict = the_table.get_celld()
    for i in range(0, len(columns)):
        # first row
        cellDict[(0, i)].set_height(.08)
        for j in range(1, len(rows) + 1):
            # values
            cellDict[(j, i)].set_height(.07)
    for i in range(1, 5):
        # first column
        cellDict[(i, -1)].set_height(.07)

    # Set text of coupled run purple
    for i in range(-1, len(columns)):
        for j in np.arange(1, 5, 2):
            cellDict[(j, i)]._text.set_color('mediumvioletred')

    # Adjust layout to make room for the table:
    plt.subplots_adjust(left=0.15, bottom=0.2)
    ax.set_title(ssp)
    plt.tight_layout()

# %% Plot Burning Embers diagramm

for ssp in Probs_const['Years'].columns.levels[0]:
    plot_burning_emeber(ssp, 2200)
    plt.tight_layout()
    plt.savefig('Plots/Comparison/5000_tip_2_5000_const_2/BurningEmber200/' + ssp + '200.png', dpi=300)

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
        if yrs_95[-1] < -10:
            print(ssp + ' ' + elem + ' ' + str(yrs_95[-1]))

    Y_diff_5[ssp] = yrs_5
    Y_diff_95[ssp] = yrs_95

print('Mean Diff [years] at 50%: ' + str(np.nanmean(list(Y_diff_5.values()))))
print('Mean Diff [years] at 95%: ' + str(np.nanmean(list(Y_diff_95.values()))))

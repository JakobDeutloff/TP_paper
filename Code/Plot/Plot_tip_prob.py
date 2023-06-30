"""
Plot Burning embers diagramm and plain tipping probabilities for any ensemble
"""
from Code.Read_data.Read_SSP_output import read_SSP_outout, read_probs
from tabulate import tabulate
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm

# %% Read data

# set ensemble
ens_name = '5000_tip'

SSP_out = read_SSP_outout(ens_name)
SSP_probs = read_probs(ens_name)


# %% Plot functions

elements = SSP_probs['Probs'].columns.levels[1].values

def plot_PvsY(ax, tip_prob, tip_prob_total, element):
    ax.plot(tip_prob['0.5'], color='k')
    ax.set_ylim(0, 1)
    ax.fill_between(tip_prob.index, tip_prob['0.05'], tip_prob['0.95'], color='grey', alpha=0.3)
    ax.fill_between(tip_prob.index, tip_prob['0.166'], tip_prob['0.833'], color='grey', alpha=0.5)
    ax.plot(tip_prob_total, color='red')
    ax.set_title(element)


def find_tip_point(tip_prob, thrsh):
    try:
        idx = np.where(tip_prob['0.5'] > thrsh)[0][0]
        year = tip_prob.index[idx]
    except:
        year = np.nan

    return year


def plot_cumm_tip(ax, tip_years, cum_tip):

    ax.set_ylim(0, 16)
    ax.plot(cum_tip)

    for element in list(tip_years.keys()):
        if not np.isnan(tip_years[element]):
            ax.text(2550, cum_tip.loc[tip_years[element]], element)

def gradientbars(bars,ydata,cmap):
    ax = bars[0].axes
    lim = ax.get_xlim()+ax.get_ylim()
    ax.axis(lim)
    for bar in bars:
        bar.set_facecolor("none")
        x, y = bar.get_xy()
        w, h = bar.get_width(), bar.get_height()
        grad = np.atleast_2d(np.linspace(0, 1*h/1, 256)).T
        #zorder of 2 to get gradients above the facecolor, but below the bar outlines
        ax.imshow(grad, extent=[x, x+w, y, y+h], origin='lower', aspect="auto", zorder=2,
                  norm=cm.colors.NoNorm(vmin=0, vmax=1), cmap=plt.get_cmap(cmap))

def plot_burning_emeber(ssp, max_year=2500):

    columns = SSP_probs['Years'].columns.levels[1]
    rows = SSP_probs['Years'].index.values
    values_max_year = SSP_out['Tip_probs_total'].loc[:max_year][ssp][columns].max().values
    values_max_total = SSP_out['Tip_probs_total'][ssp][columns].max().values

    fig, ax = plt.subplots(figsize=(8, 5))

    # set grid
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.05, 0.166, 0.5, 0.833, 0.95, 1])
    ax.grid(which='major', axis='y', linestyle='--', color='gray', zorder=0)

    # plot bars
    my_bar = ax.bar(columns, values_max_total, edgecolor='gray', zorder=3)
    # second bars to mark probs at earlier stage
    if max_year < 2500:
        # ax.bar(columns, values_max_year, edgecolor='green', color='none', zorder=3, linewidth=2)
        width = my_bar.get_children()[1].get_width()
        for i in range(len(my_bar.get_children())):
            x = my_bar.get_children()[i].get_x()
            y = values_max_year[i]
            ax.plot([x, x+width], [y, y], color='green', linewidth=2, zorder=3)

    gradientbars(my_bar, values_max_total, 'YlOrRd')
    plt.xticks([])
    first, last = ax.get_xlim()
    plt.xlim(left=first + 0.7, right=last - 0.7)
    ax.set_ylabel('Probability of Tipping')

    # draw table
    # modify table text
    text = SSP_probs['Years'][ssp][columns].values
    # Exclude values beyond max_year
    #text[text > max_year] = np.nan

    text_str = [[None for _ in range(text.shape[1])] for _ in range(text.shape[0])]
    for i in range(text.shape[0]):
        for j in range(text.shape[1]):
            try:
                text_str[i][j] = str(int(text[i, j]))
            except:
                text_str[i][j] = ''

    the_table = ax.table(cellText=text_str, rowLabels=rows, colLabels=columns,
                         loc='bottom', cellLoc='center', colLoc='center', rowLoc='center', edges='horizontal')

    # Set cell heights
    cellDict = the_table.get_celld()
    for i in range(0, len(columns)):
        # first row
        cellDict[(0, i)].set_height(.08)
        for j in range(1, len(rows) + 1):
            # values
            cellDict[(j, i)].set_height(.07)
    for i in range(1, 6):
        # first column
        cellDict[(i, -1)].set_height(.07)

    # Set text color if crossing happens after max_year
    if max_year < 2500:
        for i in range(0, len(columns)):
            for j in range(1, len(rows) + 1):
                # values
                if cellDict[(j, i)].get_text().get_text() < str(max_year):
                    cellDict[(j, i)]._text.set_color('green')


    # Adjust layout to make room for the table:
    plt.subplots_adjust(left=0.15, bottom=0.3)
    ax.set_title(ssp)

# %% Plot Burning Embers diagramm

for ssp in SSP_probs['Years'].columns.levels[0]:
    plot_burning_emeber(ssp, 2500)
    plt.savefig('Plots/SSP_ensembles/' + ens_name + '/BurningEmber/'+ssp+'.png')


# %% Plot plain tipping probabilities

for ssp in SSP_probs['Probs'].columns.levels[0].values:

    fig, axes = plt.subplots(4, 4, figsize=(10, 10), sharex='col', sharey='row')

    for i in range(len(SSP_probs['Probs'].columns.levels[1])):
        plot_PvsY(axes.flatten()[i], SSP_probs['Probs'][ssp, elements[i]], SSP_out['Tip_probs_total'][ssp, elements[i]],
                  elements[i])

    for ax in axes[:, 0]:
        ax.set_ylabel('P Tip')

    for ax in axes[3, :]:
        ax.set_xlabel('Year')

    fig.suptitle(ssp)
    plt.tight_layout()
    plt.savefig('Plots/SSP_ensembles/' + ens_name + '/P_tip/' + ssp + '.png')

# %% Plot tipping probability of only one element
ssp = 'ssp245'
element = 'AMAZ'
fig, ax = plt.subplots()
plot_PvsY(ax, SSP_probs['Probs'][ssp, element], SSP_out['Tip_probs_total'][ssp, element], element)
ax.set_ylabel('P Tip')
ax.set_xlabel('Year')
plt.savefig('Plots/SSP_enembles/10000_const/P_tip/' + ssp + element + '.png')



"""
Plot the cumulative number of TEs crossing a certain percentile of tipping probability (Fig. 5)
"""

from Code.Read_data.Read_SSP_output import read_probs
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import matplotlib.path as mpath

mpl.use('Qt5Agg')

# %%  load data
P_tip = read_probs(ens_name='coupled_ensemble')
P_tip_const = read_probs('uncoupled_ensemble')


# %% def functions
def find_y(years, year, tol, dy):
    if np.isnan(year):
        return 0

    count = 0
    for y in np.arange(year - tol, year + tol + 1, 1):
        count += years.count(y)

    return dy * count


# %% Plot

# Set plot parameters
perc = 0.5

years_tip = P_tip['Years'].loc[perc]
ssps = ['ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585']
elements = ['GRIS', 'WAIS', 'EASB', 'AWSI', 'EAIS', 'BARI', 'GLCR',
            'REEF', 'SAHL', 'BORF', 'TUND',
            'LABC', 'AMOC',
            'PFTP', 'PFAT', 'AMAZ']

# Set colors and markers
colors = {'Cry': 'blue', 'Bio': 'green', 'AOC': 'orange', 'Carbon': 'magenta'}
TE_categories = {'Cry': ['GRIS', 'WAIS', 'EASB', 'AWSI', 'EAIS', 'BARI', 'GLCR'],
                 'Bio': ['REEF', 'SAHL', 'BORF', 'TUND'],
                 'AOC': ['LABC', 'AMOC'],
                 'Carbon': ['PFTP', 'PFAT', 'AMAZ']}
# Additional marker
star = mpath.Path.unit_regular_star(6)
circle = mpath.Path.unit_circle()
cut_star = mpath.Path(
    vertices=np.concatenate([circle.vertices, star.vertices[::-1, ...]]),
    codes=np.concatenate([circle.codes, star.codes]))

# Set marker
markers = {'GRIS': '8', 'WAIS': '<', 'EASB': 'p', 'AWSI': 'v', 'EAIS': '*', 'BARI': 'D', 'GLCR': 'P',
           'REEF': '>', 'SAHL': 'H', 'BORF': 'X', 'TUND': 'd',
           'LABC': star, 'AMOC': 'h',
           'PFTP': 's', 'PFAT': '^', 'AMAZ': 'o'}

colors_ssp = {'ssp119': '#03a9d0', 'ssp126': '#193764', 'ssp245': '#f79521', 'ssp370': '#e01e27', 'ssp434': '#f0ae63',
              'ssp460': '#c28b6d', 'ssp585': '#981b1e', 'ssp534-over': '#981b1e', 'history': 'grey'}

# Set y and difference dy between two elements tipping within year +-tol
y = np.linspace(0, 1, len(ssps))
tol = 5
dy = (y[1] - y[0]) / 8

# Plot

fig, axes = plt.subplots(2, 2, figsize=(8, 8), sharex='col', sharey='row',
                         gridspec_kw={'width_ratios': [3, 1], 'wspace': 0.0, 'hspace': 0.1})

#  Plot Cum TEs --------------------------------------------------------------------------------------------------------
edges = P_tip['cumTEs'].index.values - 0.5
edges = np.append(edges, edges[-1] + 1)
for ssp in ssps:
    for ax in axes[0, :]:
        ax.stairs(P_tip['cumTEs'].loc[:, (ssp, str(perc))], edges=edges, color=colors_ssp[ssp])
        ax.stairs(P_tip_const['cumTEs'].loc[:, (ssp, str(perc))], edges=edges, color=colors_ssp[ssp], linestyle='--')

axes[0, 0].set_ylabel('Cumulative TEs with P > ' + str(perc))

#  Plot Element Symbols ------------------------------------------------------------------------------------------------
for i in range(len(ssps)):
    years = []
    for elem in elements:

        # get color
        for cat in ['Cry', 'Bio', 'AOC', 'Carbon']:
            if elem in TE_categories[cat]:
                color = colors[cat]

        # calculate y position of marker
        year = years_tip[ssps[i]][elem]
        dif_y = find_y(years, year, tol, dy)
        years.append(year)

        # write labels only for ssp585
        if ssps[i] == 'ssp585':
            for ax in axes[1, :]:
                ax.plot(year, y[i] + dif_y, marker=markers[elem], color=color, alpha=0.5, label=elem, linestyle='')
        else:
            if not perc == 0.5:  # Overlap of AMAZ for 0.5
                for ax in axes[1, :]:
                    ax.plot(year, y[i] + dif_y, marker=markers[elem], color=color, alpha=0.5, linestyle='')
            else:
                axes[1, 0].plot(year, y[i] + dif_y, marker=markers[elem], color=color, alpha=0.5, linestyle='')
                if elem == 'AMAZ' and ssps[i] == 'ssp245':
                    pass
                else: axes[1, 1].plot(year, y[i] + dif_y, marker=markers[elem], color=color, alpha=0.5, linestyle='')

# Set y ticks and grid of element plots
yticks = y - 2 * dy
yticks = np.append(yticks, y[-1] + 6 * dy)
for ax in axes[1, :]:
    ax.set_yticks(yticks, minor=True)
    ax.yaxis.grid(True, which='minor', color='grey', linestyle='--')
    ax.set_ylim([yticks[0], yticks[-1]])
    ax.set_yticks(y + 2 * dy, minor=False)
    ax.set_yticklabels(ssps, rotation=0)
    yticks_plt = ax.yaxis.get_major_ticks()
    for tck in yticks_plt:
        tck.label._color = colors_ssp[tck.label._text]

# Set xticks of plots
if perc < 0.5:
    axes[1, 0].set_xlim([1990, 2200])
    axes[1, 0].set_xticks([1990, 2050, 2100, 2150, 2200], minor=False)
else:
    axes[1, 0].set_xlim([2010, 2200])
    axes[1, 0].set_xticks([2010, 2050, 2100, 2150, 2200], minor=False)
axes[1, 1].set_xlim([2200, 2500])
axes[1, 1].set_xticks([2300, 2400, 2500], minor=False)
fig.text(0.52, 0.125, 'Year')

# hide the spines between horizontal axes
for ax in axes[:, 0]:
    ax.spines.right.set_visible(False)
    ax.yaxis.tick_left()

for ax in axes[:, 1]:
    ax.spines.left.set_visible(False)
    ax.yaxis.tick_right()
    ax.tick_params(labelleft=False)

# Vertical lines
for ax in axes.flatten():
    for year in [2050, 2100, 2200]:
        ax.axvline(year, linestyle='--', linewidth=0.7, color='grey')

# Set legend
fig2, ax = plt.subplots()  # Fake plot for add. handles and labels
ax.plot([1, 2, 3], [2, 3, 4], linestyle='--', color='k', label='Uncoupled')
ax.plot([1, 2, 3], [2, 3, 4], linestyle='-', color='k', label='Coupled')
hd, lb = ax.get_legend_handles_labels()
fig2.clf()
ax.cla()

handles, labels = axes[1, 0].get_legend_handles_labels()
handles.append(hd[0])
labels.append(lb[0])
handles.append(hd[1])
labels.append(lb[1])

fig.subplots_adjust(bottom=0.17)
fig.legend(handles, labels, bbox_to_anchor=(0.93, 0.1), ncol=6)
fig.savefig('Plots/Comparison/5000_tip_2_5000_const_2/Cum_TEs/perc' + str(perc) + '_cum_TEs.png', dpi=300,
            bbox_inches='tight')

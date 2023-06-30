"""
Plot results from variance analysis comparing LH and random sampling
"""
import matplotlib.pyplot as plt
import numpy as np
from Code.Read_data.read_var_and_stab_analysis import LH585, Rand585

# %%
fig, axes = plt.subplots(1, 3, figsize=(8, 4), sharex='col')

elements = ['PFTP', 'AMAZ', 'PFAT']


def plot_var(ax, elem):
    # calculate mean of individual ensembles
    mean_LH = LH585.xs((elem, 'C_stock'), level=(2, 3), axis=1).groupby(level=0, axis=1).mean()
    mean_Rand = Rand585.xs((elem, 'C_stock'), level=(2, 3), axis=1).groupby(level=0, axis=1).mean()
    # calculate mean over ensembles and std
    ens_mean_LH = mean_LH.mean(axis=1)
    ens_std_LH = mean_LH.std(axis=1)
    ens_mean_Rand = mean_Rand.mean(axis=1)
    ens_std_Rand = mean_Rand.std(axis=1)

    # Plot mean and std
    ax.plot(ens_mean_Rand, color='blue', label=r'$\mu$ Random', linestyle='-')
    ax.fill_between(ens_mean_Rand.index, ens_mean_Rand - 50 * ens_std_Rand, ens_mean_Rand + 50 * ens_std_Rand,
                    color='blue',
                    alpha=0.2, label=r'$\mathrm{50 \cdot \sigma ~ Random}$')
    ax.plot(ens_mean_LH, color='red', label=r'$\mu$ LH', linestyle='--')
    ax.fill_between(ens_mean_LH.index, ens_mean_LH - 50 * ens_std_LH, ens_mean_LH + 50 * ens_std_LH, color='red',
                    alpha=0.5, label=r'$\mathrm{50 \cdot \sigma ~ LH}$')

    ax.set_title(elem)
    ax.set_xlabel('Year')


plot_var(axes[0], 'PFTP')
plot_var(axes[1], 'AMAZ')
plot_var(axes[2], 'PFAT')
axes[0].set_ylabel('Cumulative C Emissions [GtC]')

# Legend
handles, labels = axes[0].get_legend_handles_labels()
fig.subplots_adjust(bottom=0.2)
fig.legend(handles, labels, ncol=4, loc='lower center')

plt.savefig('Plots/Calibration/LHvsRandom.png', bbox_inches='tight', dpi=500)
plt.show()

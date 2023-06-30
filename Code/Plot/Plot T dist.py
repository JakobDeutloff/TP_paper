"""
Script to plot histograms of T and calculate quantiles and probabilities of exceedance.
"""
from Code.Read_data.Read_SSP_output import read_res_years, read_T_ens
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from scipy.stats import t
mpl.use('Qt5Agg')


# %% define plot function

def plot_hist(ssp, year, ax):
    val_ens = {}
    val_tip = res_tip[ssp]['T'].loc[year].values
    for i in range(9):
        val = T_ens[i][ssp]['T'].loc[year].values
        val_ens[str(i)] = val
    val_ens['9'] = val_tip

    val_const = res_const[ssp]['T'].loc[year].values

    val_tip_unpacked = np.concatenate(list(val_ens.values()))
    max = np.ceil(np.max(val_tip_unpacked))
    min = 0
    step = 0.1
    bins = np.arange(min, max + step, step)

    # Plot hist
    hist_tip = ax.hist(val_tip_unpacked, density=True, bins=bins, color='red', alpha=0.5, zorder=2)
    ax.hist(val_const, density=True, bins=bins, color='k', alpha=0.5, zorder=1)

    # Plot tail
    axins = ax.inset_axes([0.5, 0.5, 0.47, 0.47])
    axins.hist(val_tip_unpacked, density=True, bins=bins, color='red', alpha=0.5, zorder=2)
    axins.hist(val_const, density=True, bins=bins, color='k', alpha=0.5, zorder=1)

    # Set Labels
    ax.set_xlabel('GMST anomaly [°C]')
    ax.set_title(ssp)

    return axins, ax, hist_tip, val_ens, val_const


# Function for calculation of exceedence probs
def calc_P(temp, val_ens, val_const):
    # coupled
    p_tip = []
    for i in range(10):
        n_tip = np.sum(val_ens[str(i)] > temp)
        p_tip += [n_tip / len(val_ens[str(i)])]

    p_tip_mean = np.mean(p_tip)
    p_tip_std = np.std(p_tip)
    print('Probability from coupled run of exceeding ' + str(temp) + '°C : ' + str(p_tip_mean) + '+-' + str(p_tip_std))

    # uncoupled
    n_const = np.sum(val_const > temp)
    p_const = n_const / len(val_const)
    print('Probability from uncoupled run of exceeding ' + str(temp) + '°C : ' + str(p_const))
    print ('P-value: ' + str(test_significance(p_tip_mean, p_const, p_tip_std, 10)))

def calc_quants(val_ens, val_const):
    # coupled
    lower = []
    median = []
    upper = []

    for i in range(10):
        lower += [np.quantile(val_ens[str(i)], 0.05)]
        median += [np.quantile(val_ens[str(i)], 0.5)]
        upper += [np.quantile(val_ens[str(i)], 0.95)]

    print('Quantiles coupled ens: 0.05:' + str(np.mean(lower).round(3)) + '+-' + str(np.std(lower).round(3)) +
          ', 0.5: ' + str(np.mean(median).round(3)) + '+-' + str(np.std(median).round(3)) +
          ', 0.95: ' + str(np.mean(upper).round(3)) + '+-' + str(np.std(upper).round(3)))

    # uncoupled
    print('Quantiles uncoupled ens: 0.05:' + str(np.quantile(val_const, 0.05).round(3)) +
          ', 0.5: ' + str(np.quantile(val_const, 0.5).round(3)) +
          ', 0.95: ' + str(np.quantile(val_const, 0.95).round(3)))

    # Test significance
    p_vals = []
    p_vals += [test_significance(np.mean(lower), np.quantile(val_const, 0.05), np.std(lower), 10)]
    p_vals += [test_significance(np.mean(median), np.quantile(val_const, 0.5), np.std(median), 10)]
    p_vals += [test_significance(np.mean(upper), np.quantile(val_const, 0.95), np.std(upper), 10)]
    print('P-value: 0.05:' + str(p_vals[0]) + ', 0.5: ' + str(p_vals[1]) + ', 0.95: ' + str(p_vals[2]))


def test_significance(mu_c, mu_u, std_c, n_c):

    t_val = np.sqrt(n_c) * (mu_c - mu_u) / std_c
    p = 1 - t.cdf(t_val, df=n_c-1)
    return p

# %% Load Data

res_tip = read_res_years('5000_tip_2')
res_const = read_res_years('5000_const_2')
T_ens = read_T_ens()

# %% plot ssp 119 and ssp 245 2100


fig, axes = plt.subplots(1, 2, figsize=(10, 5))
year = 2100

# SSP119
axins, ax, hist_tip, val_ens_119, val_const_119 = plot_hist('ssp119', year, axes[0])

# get and set limits
x_max = 3.5
x_min = 2
y_min = 0
idx_ymax = (np.abs(hist_tip[1] - x_min)).argmin()
y_max = hist_tip[0][idx_ymax]
axins.set_xlim(x_min, x_max)
axins.set_ylim(y_min, y_max)
axes[0].set_ylabel('Normalized Frequency')
ax.set_xlim(0, 5)

# SSP 245
axins, ax, hist_tip, val_ens_245, val_const_245 = plot_hist('ssp245', year, axes[1])

# get and set limits
x_max = 6
x_min = 4
y_min = 0
idx_ymax = (np.abs(hist_tip[1] - x_min)).argmin()
y_max = hist_tip[0][idx_ymax]
axins.set_xlim(x_min, x_max)
axins.set_ylim(y_min, y_max)
ax.set_xlim(0, 14)

fig.tight_layout()
fig.savefig('Plots/Comparison/5000_tip_2_5000_const_2/T_hists/ssp119_245' + str(year) + 'hist.png', dpi=300)
fig.show()

#  Calculate probabilies of exeedence (mean +- one std)
calc_P(2, val_ens_119, val_const_119)
calc_P(4, val_ens_245, val_const_245)

# %% Calculate Quantiles
calc_quants(val_ens_119, val_const_119)
calc_quants(val_ens_245, val_const_245)

# %% plot ssp 119 and ssp 245 2300


fig, axes = plt.subplots(1, 2, figsize=(10, 5))
year = 2300

# SSP119
axins, ax, hist_tip, val_ens_119, val_const_119 = plot_hist('ssp119', year, axes[0])

# get and set limits
x_max = 3.5
x_min = 2
y_min = 0
idx_ymax = (np.abs(hist_tip[1] - x_min)).argmin()
y_max = hist_tip[0][idx_ymax]
axins.set_xlim(x_min, x_max)
axins.set_ylim(y_min, y_max)
axes[0].set_ylabel('Normalized Frequency')
ax.set_xlim(0, 5)

# SSP 245
axins, ax, hist_tip, val_ens_245, val_const_245 = plot_hist('ssp245', year, axes[1])

# get and set limits
x_max = 8
x_min = 5
y_min = 0
idx_ymax = (np.abs(hist_tip[1] - x_min)).argmin()
y_max = hist_tip[0][idx_ymax]
axins.set_xlim(x_min, x_max)
axins.set_ylim(y_min, y_max)
ax.set_xlim(0, 15)

fig.tight_layout()
fig.savefig('Plots/Comparison/5000_tip_2_5000_const_2/T_hists/ssp119_245' + str(year) + 'hist.png', dpi=300)
fig.show()

#  Calculate probabilies of exeedence (mean +- one std)
calc_P(2, val_ens_119, val_const_119)
calc_P(4, val_ens_245, val_const_245)

# %% Calculate Quantiles
calc_quants(val_ens_119, val_const_119)
calc_quants(val_ens_245, val_const_245)

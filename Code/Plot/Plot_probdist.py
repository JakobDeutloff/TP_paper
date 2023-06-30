"""
Plots the probability distributions used in the analysis (Fig. 2, Fig. S5-S7)
"""
import matplotlib.pyplot as plt
import numpy as np
from Code.Calibrate.find_distributions import P_log_norm_dist, P_perf_dists, Ts_perf_dists, K_beta_dists, F_beta_dists
from Code.Read_data.Constants_new import P_min, P_max, P_mean, par_P_perf_dist, Ts_max, Ts_min, Ts_mean, r_PFTP, K_min, \
    K_max, F_min, F_max, par_Ts_perf_dist

elements = ['PFTP', 'AMOC', 'GRIS', 'WAIS', 'AMAZ', 'BORF', 'AWSI', 'EAIS', 'EASB', 'GLCR', 'LABC', 'TUND', 'PFAT',
            'BARI', 'REEF', 'SAHL']
ts_elements = ['PFTP', 'AMAZ']

# %% calculate RMSE
RMSE_P = {}
RMSE_Ts = {}

for elem in elements:
    RMSE_P[elem] = np.sqrt(((P_perf_dists[elem].ppf(0.05) - P_min[elem]) ** 2 + \
                            (P_perf_dists[elem].ppf(0.5) - P_mean[elem]) ** 2 + \
                            (P_perf_dists[elem].ppf(0.95) - P_max[elem]) ** 2) / 3)
    if isinstance(RMSE_P[elem], type(np.array([0]))):
        RMSE_P[elem] = RMSE_P[elem][0]

for elem in ts_elements:
    RMSE_Ts[elem] = np.sqrt(((Ts_perf_dists[elem].ppf(0.05) - Ts_min[elem]) ** 2 + \
                             (Ts_perf_dists[elem].ppf(0.5) - Ts_mean[elem]) ** 2 + \
                             (Ts_perf_dists[elem].ppf(0.95) - Ts_max[elem]) ** 2) / 3)
    if isinstance(RMSE_P[elem], type(np.array([0]))):
        RMSE_Ts[elem] = RMSE_Ts[elem][0]


# %% Plot H of PFTP and AMAZ distributions

def plot_ts_dist(ax, elem):
    x = np.arange(0, 351)
    y = Ts_perf_dists[elem].cdf(x=x)
    ax.plot(x, y, color='blue')
    ax.plot([Ts_min[elem], Ts_mean[elem], Ts_max[elem]], [0.05, 0.5, 0.95], color='r', linestyle='', marker='o',
            markersize=5)
    ax.set_title(elem + ', ' + par_Ts_perf_dist['PFTP']['name'] + ', RMSE=' + str(RMSE_Ts[elem].round(2)))
    ax.set_xlabel('Timescale [Years]')


fig, axes = plt.subplots(1, 2, sharey='row', figsize=(8, 2.5))

for i in range(len(Ts_perf_dists)):
    plot_ts_dist(axes[i], ts_elements[i])

axes[0].set_ylabel('Cumulative Probability')
plt.tight_layout()
plt.savefig('Plots/Probability_dists/Ts_prob_dist.png', dpi=500)
plt.show()


# %% Plot P for all TEs
def plot_distr_perf(elem, ax):

    x = np.linspace(0, 10, 200)
    y = P_perf_dists[elem].cdf(x=x)
    ax.plot(x, y, color='blue')
    ax.plot([P_min[elem], P_mean[elem], P_max[elem]], [0.05, 0.5, 0.95], color='r', linestyle='', marker='o')
    ax.set_title(elem + ', ' + par_P_perf_dist[elem]['name'] + ', RMSE=' + str(np.round(RMSE_P[elem], decimals=2)))


fig, axes = plt.subplots(4, 4, sharex='col', sharey='row', figsize=(11, 9))
for i in range(len(elements)):
    plot_distr_perf(elements[i], axes.flat[i])

for ax in axes[-1, :]:
    ax.set_xlabel('GMST anomaly [°C]')

for ax in axes[:, 0]:
    ax.set_ylabel('Cumulative Probability')

plt.tight_layout()
plt.savefig('Plots/Probability_dists/all_perf.png', dpi=500)
plt.show()

# %% Plot K for carbon TEs
fig, axes = plt.subplots(1, 3, figsize=(8, 3), sharey='row')

# PFTP
x1 = np.linspace(100, 275, 1000)
y1 = K_beta_dists['PFTP'].cdf(x1)
axes[0].plot(x1, y1)
axes[0].set_title('PFTP')
axes[0].axvline(K_min['PFTP'], alpha=1, color='r', linestyle='--')
axes[0].axvline(K_max['PFTP'], alpha=1, color='r', linestyle='--')
# AMAZ
x1 = np.linspace(20, 85, 10000)
y1 = K_beta_dists['AMAZ'].cdf(x1)
axes[1].plot(x1, y1)
axes[1].set_title('AMAZ')
axes[1].axvline(K_min['AMAZ'], alpha=1, color='r', linestyle='--')
axes[1].axvline(K_max['AMAZ'], alpha=1, color='r', linestyle='--')
# PFAT
x1 = np.linspace(45, 215, 10000)
y1 = K_beta_dists['PFAT'].cdf(x1)
axes[2].plot(x1, y1)
axes[2].set_title('PFAT')
axes[2].axvline(K_min['PFAT'], alpha=1, color='r', linestyle='--')
axes[2].axvline(K_max['PFAT'], alpha=1, color='r', linestyle='--')

for ax in axes:
    ax.set_xlabel('Impact [GtC]')
axes[0].set_ylabel('Cumulative Probability')

plt.tight_layout()
plt.savefig('Plots/Probability_dists/K_prob_dist.png', dpi=500)

plt.show()

# %% Plot F for PFAT

fig, axes = plt.subplots(1, 2, figsize=(8, 3), sharey='row')

# PFAT100
x = np.linspace(4, 16, 10000)
y = F_beta_dists['PFAT100'].cdf(x)
axes[0].plot(x, y)
axes[0].set_title(r'$F_{\mathrm{100}}$')
axes[0].axvline(F_min['PFAT100'], alpha=1, color='r', linestyle='--')
axes[0].axvline(F_max['PFAT100'], alpha=1, color='r', linestyle='--')

# PFAT300
x = np.linspace(10, 40, 10000)
y = F_beta_dists['PFAT300'].cdf(x)
axes[1].plot(x, y)
axes[1].set_title(r'$F_{\mathrm{300}}$')
axes[1].axvline(F_min['PFAT300'], alpha=1, color='r', linestyle='--')
axes[1].axvline(F_max['PFAT300'], alpha=1, color='r', linestyle='--')

for ax in axes:
    ax.set_xlabel('Feedback [GtC/°C]')
axes[0].set_ylabel('Cumulative Probability')

plt.tight_layout()
plt.savefig('Plots/Probability_dists/F_prob_dist.png', dpi=500)
plt.show()


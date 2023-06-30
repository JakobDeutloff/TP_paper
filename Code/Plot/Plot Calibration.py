"""
Plot the calibrated r of PFTP and AMAZ (Fig. S2) as well as S for AMAZ (Fig. 1) under SSP5-8.5
"""

import pandas as pd
import matplotlib.pyplot as plt
from Code.Read_data.Constants_new import Ts_min, Ts_max, Ts_mean
import pickle
from Code.Calibrate.Functions import find_threshold, find_trg_year

# %% Read data
r_PFTP = pd.read_csv('Data/Calibration/R_Ts_E_PFTP.csv', index_col=[0])
r_AMAZ = pd.read_csv('Data/Calibration/R_Ts_E_AMAZ.csv', index_col=[0])
mean_AMAZ = pickle.load(open('Data/Calibration/mean_AMAZ.pkl', 'rb'))

# %% Plot calibrated r

fig, axes = plt.subplots(1, 2, figsize=(8, 2.5))

Ts = [Ts_min, Ts_mean, Ts_max]

axes[0].plot(r_PFTP.timescale, r_PFTP.rate, color='k')
axes[0].set_title('PFTP')
axes[0].set_ylabel(r'Rate [$\mathrm{yr^{-1}}$]')

axes[1].plot(r_AMAZ.timescale, r_AMAZ.rate, color='k')
axes[1].set_title('AMAZ')


for t in Ts:
    axes[0].axvline(t['PFTP'], color='grey', linestyle='--')
    axes[1].axvline(t['AMAZ'], color='grey', linestyle='--')

for ax in axes:
    ax.set_xlabel('Timescale [Years]')
plt.tight_layout()

plt.savefig('Plots/Calibration/Ts_PFTP_AMAZ.png', dpi=500)
plt.show()

# %% plot S of AMAZ

trg_year = find_trg_year(mean_AMAZ['tipping_elements'].xs(('AMAZ', 'bool_tip'), level=(1, 2), axis=1))
thr_year = find_threshold(mean_AMAZ['tipping_elements'].xs(('AMAZ', 'C_stock'), level=(1, 2), axis=1), 0.995)

fig, ax = plt.subplots(figsize=(4, 2.5))
mean_AMAZ['tipping_elements'].loc[2050:2200, ('ssp585', 'AMAZ', 'C_stock')].plot(ax=ax, color='k')
ax.axvspan(trg_year, thr_year, alpha=0.2, color='red')
ax.set_ylabel(r'$S$ [GtC]')
ax.set_xlabel('Year')
plt.tight_layout()

plt.savefig('Plots/Calibration/mean_Ts_AMAZ.png', dpi=500)
plt.show()

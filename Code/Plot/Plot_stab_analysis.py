"""
Plot stability analysis results for CTEM output and GMST
"""
from Code.Read_data.read_var_and_stab_analysis import T
#%%
import matplotlib.pyplot as plt
import numpy as np

# %% Plot temperature
fig, ax = plt.subplots(figsize=(5, 4))
mean = T.xs('0.5', level=1, axis=1).mean(axis=1)
std = T.xs('0.5', level=1, axis=1).std(axis=1)
mean05 = T.xs('0.05', level=1, axis=1).mean(axis=1)
std05 = T.xs('0.05', level=1, axis=1).std(axis=1)
mean95 = T.xs('0.95', level=1, axis=1).mean(axis=1)
std95 = T.xs('0.95', level=1, axis=1).std(axis=1)

ax.fill_between(mean05.index, mean05 - 2*std05, mean05 + 2*std05,
                color='k', alpha=0.5)
ax.plot(mean05, color='k', linewidth=0.7, linestyle='--', label='5th percentile')

ax.fill_between(mean95.index, mean95 - 2*std95, mean95 + 2*std95,
                color='k', alpha=0.5)
ax.plot(mean95, color='k', linewidth=0.7, label='95th percentile')

ax.fill_between(mean.index, mean - 2*std, mean + 2*std,
                color='k', alpha=0.5, label=r'$\mathrm{2 \cdot \sigma} $')
ax.plot(mean, color='r', linewidth=0.7, label='Median')


ax.set_xlabel('Year')
ax.set_ylabel('GMST anomaly [°C]')
ax.legend()
plt.tight_layout()
plt.savefig('Plots/Calibration/stab_analysis/T_stab.png', dpi=200)
plt.show()


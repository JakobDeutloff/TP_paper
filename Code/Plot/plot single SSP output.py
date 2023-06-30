"""
Old script to plot T and C from SSP ensembles - needs updating
"""
import matplotlib.pyplot as plt
from Code.Read_data.Read_SSP_output import read_SSP_outout
import numpy as np

# %% choose ssps and colors
colors = {'ssp119': '#03a9d0', 'ssp126': '#193764', 'ssp245': '#f79521', 'ssp370': '#e01e27', 'ssp434': '#f0ae63',
          'ssp460': '#c28b6d', 'ssp585': '#981b1e', 'ssp534-over': '#981b1e', 'history': 'grey'}
tier1_ssps = ['ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585']
tier2_ssps = ['ssp434', 'ssp460', 'ssp534-over']

# %% Import Data
ens_name = '5000_const'
res = read_SSP_outout(ens_name)
T_const = res['T']
C_const = res['C']
# %% Define ssp plot function

def plot_T_ssp(ssp, ax, quants, multiple=True):

    ax.plot(T_const.loc[2014:2500, (ssp, '0.5')], color=colors[ssp], label=ssp)

    if '0.95' in quants:
        ax.fill_between(np.arange(2016, 2501), *T_const.loc[2016:2500, (ssp, ['0.05', '0.95'])].values.T, color=colors[ssp],
                        alpha=0.2, lw=0)
        ax.plot(np.arange(2016, 2501), T_const.loc[2016:2500, (ssp, ['0.05', '0.95'])].values, color=colors[ssp], alpha=0.5,
                lw=1)
    if '0.833' in quants:
        ax.fill_between(np.arange(2016, 2501), *T_const.loc[2016:2500, (ssp, ['0.166', '0.833'])].values.T, color=colors[ssp],
                        alpha=0.2, lw=0)
        ax.plot(np.arange(2016, 2501), T_const.loc[2016:2500, (ssp, ['0.166', '0.833'])].values, color=colors[ssp], alpha=0.5,
                lw=1)

    # historical
    if multiple:
        ax.plot(T_const.loc[1850:2016, ('ssp245', '0.5')], color=colors['history'], label='historical')
        ax.fill_between(np.arange(1850, 2017), *T_const.loc[1850:2016, ('ssp245', ['0.05', '0.95'])].values.T,
                        color=colors['history'], alpha=0.2, lw=0)
        ax.plot(np.arange(1850, 2017), T_const.loc[1850:2016, ('ssp245', ['0.05', '0.95'])].values, color=colors['history'],
                alpha=0.5, lw=1)
    ax.set_title(ssp)

# Define C concentration plot function
def plot_C_ssp(ssp, ax, quants, multiple=True):

    ax.plot(C_const.loc[2016:2500, (ssp, 'carbon_dioxide', '0.5')], color='green', label='CO2')

    if '0.95' in quants:
        ax.fill_between(np.arange(2016, 2501), *C_const.loc[2016:2500, (ssp, 'carbon_dioxide', ['0.05', '0.95'])].values.T, color='green',
                        alpha=0.2, lw=0)
        ax.plot(np.arange(2016, 2501), C_const.loc[2016:2500, (ssp, 'carbon_dioxide', ['0.05', '0.95'])].values, color='green', alpha=0.5,
                lw=1)
    if '0.833' in quants:
        ax.fill_between(np.arange(2016, 2501), *C_const.loc[2016:2500, (ssp, 'carbon_dioxide', ['0.166', '0.833'])].values.T, color='green',
                        alpha=0.2, lw=0)
        ax.plot(np.arange(2016, 2501), C_const.loc[2016:2500, (ssp, 'carbon_dioxide', ['0.166', '0.833'])].values, color='green', alpha=0.5,
                lw=1)
    # historical
    if multiple:
        ax.plot(C_const.loc[1850:2016,  ('ssp245', 'carbon_dioxide', '0.5')], color=colors['history'], label='historical')
        ax.fill_between(np.arange(1850, 2017), *C_const.loc[1850:2016, ('ssp245', 'carbon_dioxide', ['0.05', '0.95'])].values.T,
                        color=colors['history'], alpha=0.2, lw=0)
        ax.plot(np.arange(1850, 2017), C_const.loc[1850:2016, ('ssp245', 'carbon_dioxide', ['0.05', '0.95'])].values, color=colors['history'],
                alpha=0.5, lw=1)
    ax.set_title(ssp)

# Def Methane conc. plot function
def plot_m_ssp(ssp, ax, quants, multiple=True):

    ax.plot(C_const.loc[2016:2500, (ssp, 'methane', '0.5')], color='darkslategray', label='Methane')

    if '0.95' in quants:
        ax.fill_between(np.arange(2016, 2501), *C_const.loc[2016:2500, (ssp, 'methane', ['0.05', '0.95'])].values.T, color='darkslategray',
                        alpha=0.2, lw=0)
        ax.plot(np.arange(2016, 2501), C_const.loc[2016:2500, (ssp, 'methane', ['0.05', '0.95'])].values, color='darkslategray', alpha=0.5,
                lw=1)
    if '0.833' in quants:
        ax.fill_between(np.arange(2016, 2501), *C_const.loc[2016:2500, (ssp, 'methane', ['0.166', '0.833'])].values.T, color='darkslategray',
                        alpha=0.2, lw=0)
        ax.plot(np.arange(2016, 2501), C_const.loc[2016:2500, (ssp, 'methane', ['0.166', '0.833'])].values, color='darkslategray', alpha=0.5,
                lw=1)
    # historical
    if multiple:
        ax.plot(C_const.loc[1850:2016,  ('ssp245', 'methane', '0.5')], color=colors['history'], label='historical')
        ax.fill_between(np.arange(1850, 2017), *C_const.loc[1850:2016, ('ssp245', 'methane', ['0.05', '0.95'])].values.T,
                        color=colors['history'], alpha=0.2, lw=0)
        ax.plot(np.arange(1850, 2017), C_const.loc[1850:2016, ('ssp245', 'methane', ['0.05', '0.95'])].values, color=colors['history'],
                alpha=0.5, lw=1)
    ax.set_title(ssp)
    ax.set_ylim(800, 3500)

# %% Plot all ssps seperately
fig, axes = plt.subplots(3, 2, figsize=(10, 12), sharex='col', sharey=True)

for i in range(len(tier1_ssps)):
    plot_T_ssp(tier1_ssps[i], axes.flatten()[i], quants=['0.95', '0.833'])

for ax in axes[:, 0]:
    ax.set_ylabel('temperature anomaly relative to 1850-1900 / K')

for ax in axes[2, :]:
    ax.set_xlabel('Year')

axes[1, 1].set_xlabel('Year')
axes[2, 1].remove()


plt.tight_layout()
plt.savefig('Plots/SSP_ensembles/' + ens_name + '/T_const_split_tier1.png')
plt.show()

# %% Plot only one ssp
ssp = 'ssp245'
fig, ax = plt.subplots()
plot_T_ssp(ssp, ax, quants=['0.95', '0.833'])
ax.set_ylabel('temperature anomaly relative to 1850-1900 / K')
ax.set_xlabel('Year')
plt.savefig('Plots/SSP_enembles/' + ssp + '.png')
plt.show()
# %% plot alls SSP_enembles in one plot 0.166 - 0.833 percentile

fig, ax = plt.subplots(figsize=(5, 4))

for i in range(len(tier1_ssps)):
    plot_T_ssp(tier1_ssps[i], ax, quants=['0.833'], multiple=False)


ax.plot(T_const.loc[1850:2016, ('ssp245', '0.5')], color=colors['history'], label='historical')
ax.fill_between(np.arange(1850, 2017), *T_const.loc[1850:2016, ('ssp245', ['0.166', '0.833'])].values.T,
                color=colors['history'], alpha=0.2, lw=0)
ax.plot(np.arange(1850, 2017), T_const.loc[1850:2016, ('ssp245', ['0.166', '0.833'])].values, color=colors['history'],
        alpha=0.5, lw=1)

ax.set_ylabel('GMST anomaly [°C]')
ax.set_xlabel('Year')
ax.set_title('')
ax.legend()
plt.tight_layout()
plt.savefig('Plots/SSP_ensembles/' + ens_name + '/T_tier1_16_83.png', dpi=500)
plt.show()

# %% plot alls SSP_enembles in one plot 0.05-0.95 percentile

fig, ax = plt.subplots(figsize=(8, 7))

for i in range(len(tier1_ssps)):
    plot_T_ssp(tier1_ssps[i], ax, quants=['0.95'], multiple=False)

ax.plot(T_const.loc[1850:2016, ('ssp245', '0.5')], color=colors['history'], label='historical')
ax.fill_between(np.arange(1850, 2017), *T_const.loc[1850:2016, ('ssp245', ['0.05', '0.95'])].values.T,
                color=colors['history'], alpha=0.2, lw=0)
ax.plot(np.arange(1850, 2017), T_const.loc[1850:2016, ('ssp245', ['0.05', '0.95'])].values, color=colors['history'],
        alpha=0.5, lw=1)

ax.set_ylabel('temperature anomaly relative to 1850-1900 / K')
ax.set_xlabel('Year')
ax.set_title('0.05 - 0.95 percentile')
ax.legend()
plt.tight_layout()
plt.savefig('Plots/SSP_ensembles/' + ens_name + '/T_05-95_tier1.png')
plt.show()

# %% plot GHG concentraions
fig, axes = plt.subplots(3, 2, figsize=(10, 12), sharex='col', sharey=True)
axes2 = np.empty((3, 2), dtype=type(axes[0, 0]))
for i in range(3):
    for j in range(2):

        axes2[i, j] = axes[i, j].twinx()

for i in range(len(ssps)):
    plot_C_ssp(ssps[i], axes.flatten()[i], quants=['0.95', '0.833'])
    plot_m_ssp(ssps[i], axes2.flatten()[i], quants=['0.95', '0.833'])

for ax in axes[:, 0]:
    ax.set_ylabel(r'Atmospheric $\mathrm{CO_{2}}$ Concentration [ppm]', color='green')

for ax in axes2[:, 1]:
    ax.set_ylabel(r'Atmospheric $\mathrm{CH_{4}}$ Concentration [ppb]', color='darkslategray')
axes2[2, 0].set_ylabel(r'Atmospheric $\mathrm{CH_{4}}$ Concentration [ppb]', color='darkslategray')

for ax in axes[2, :]:
    ax.set_xlabel('Year')
axes[1, 1].set_xlabel('Year')

axes[2, 1].remove()
axes2[2, 1].remove()

plt.tight_layout()
plt.savefig('Plots/SSP_enembles/C_const_split.png')
plt.show()


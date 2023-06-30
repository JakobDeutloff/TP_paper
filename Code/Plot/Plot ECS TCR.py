import matplotlib.pyplot as plt
import numpy as np

# %% Estimates of ECS and TCR

ECS_IPCC = {'0.05': 2, '0.166': 2.5, '0.5': 3, '0.833': 4, '0.95': 5}
TCR_IPCC = {'0.05': 1.2, '0.166': 1.4, '0.5': 1.8, '0.833': 2.2, '0.95': 2.4}

ECS_CMIP6 = np.array(
    [4.7, 3.9, 3.2, 3, 3.3, 2.3, 5.2, 4.8, 4.8, 4.3, 4.8, 5.6, 5.3, 4.3, 4.3, 3, 3.9, 2.6, 2.7, 3.1, 2.4, 5.6,
     5.4, 1.8, 1.9, 4.6, 4.5, 3.7, 2.7, 2.6, 3.0, 3.0, 3.2, 4.7, 2.5, 3.7, 5.3])
TCR_CMIP6 = np.array(
    [2.1, 2, 2, 1.7, 1.8, 1.7, 2, 2, 2.1, 2.5, 1.9, 2.7, 3, 2.6, 2.1, 2.1, 1.6, 1.8, 1.9, 1.7, 2.6, 2.6, 1.7,
     1.3, 2.3, 1.4, 1.9, 1.6, 1.6, 1.7, 1.8, 1.6, 2.7, 1.6, 1.5, 2.3, 2.8])
# Calculate percentiles of CMIP6
TCR_CMIP6_q = {'0.05': np.quantile(TCR_CMIP6, 0.05), '0.166': np.quantile(TCR_CMIP6, 0.166),
               '0.5': np.quantile(TCR_CMIP6, 0.5), '0.833': np.quantile(TCR_CMIP6, 0.833),
               '0.95': np.quantile(TCR_CMIP6, 0.95)}
ECS_CMIP6_q = {'0.05': np.quantile(ECS_CMIP6, 0.05), '0.166': np.quantile(ECS_CMIP6, 0.166),
               '0.5': np.quantile(ECS_CMIP6, 0.5), '0.833': np.quantile(ECS_CMIP6, 0.833),
               '0.95': np.quantile(ECS_CMIP6, 0.95)}

ECS_FaIR = {'0.05': 1.94, '0.166': 2.36, '0.5': 3.24, '0.833': 4.74, '0.95': 6.59}
TCR_FaIR = {'0.05': 1.3, '0.166': 1.48, '0.5': 1.79, '0.833': 2.15, '0.95': 2.44}

ECS = {'CMIP6': ECS_CMIP6_q, 'IPCC AR6': ECS_IPCC, 'FaIR': ECS_FaIR}
TCR = {'CMIP6': TCR_CMIP6_q, 'IPCC AR6': TCR_IPCC, 'FaIR': TCR_FaIR}
# %% Plot estimates

fig, axes = plt.subplots(1, 2, figsize=(10, 5))

names = list(ECS.keys())
x = [1, 2, 3]
width = 0.5

for i in range(3):

    # ECS
    axes[0].bar(x[i], height=ECS[names[i]]['0.95']-ECS[names[i]]['0.05'], width=width, bottom=ECS[names[i]]['0.05'],
                color='white', edgecolor='k')
    axes[0].bar(x[i], height=ECS[names[i]]['0.833']-ECS[names[i]]['0.166'], width=width, bottom=ECS[names[i]]['0.166'],
                color='grey', alpha=0.8)
    axes[0].plot([x[i]-(width/2)+0.01, x[i]+(width/2)-0.01], [ECS[names[i]]['0.5'], ECS[names[i]]['0.5']],
                 color='k', label='Median')

    # TCR
    axes[1].bar(x[i], height=TCR[names[i]]['0.95']-TCR[names[i]]['0.05'], width=width, bottom=TCR[names[i]]['0.05'],
                color='white', edgecolor='k', label='Very Likely Range')
    axes[1].bar(x[i], height=TCR[names[i]]['0.833']-TCR[names[i]]['0.166'], width=width, bottom=TCR[names[i]]['0.166'],
                color='grey', alpha=0.8, label='Likely Range')
    axes[1].plot([x[i]-(width/2)+0.01, x[i]+(width/2)-0.01], [TCR[names[i]]['0.5'], TCR[names[i]]['0.5']],
                 color='k', label='Median')

for ax in axes:
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(names)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
axes[0].set_ylabel('ECS [°C]')
axes[0].set_ylim([1.5, 7])
axes[1].set_ylabel('TCR [°C]')
axes[1].set_ylim([1.1, 2.8])

handles, labels = axes[1].get_legend_handles_labels()
fig.legend(handles[2:5], labels[2:5], ncol=3, bbox_to_anchor=(0.75, -0.001))

fig.savefig('Plots/ECS_TCR.png', dpi=300, bbox_inches='tight')

plt.show()

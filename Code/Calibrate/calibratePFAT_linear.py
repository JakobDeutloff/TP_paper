"""
Script to check implementation and plot PFAT
"""

from Code.Read_data.read_simplified_RCMIP import *
from Code.FAIR2.fair.fair_runnter_intanal import run_FaIR_intanal
import matplotlib.pyplot as plt
from Code.Read_data.Constants_new import *


# %% get emissions
ssp = 'ssp460'
emms = SSP_emms[[ssp]]
emms.columns = emms.columns.remove_unused_levels()
forc = SSP_forc[[ssp]]
forc.columns = forc.columns.remove_unused_levels()

# %% get sample of tipping probs, timescales and impacts of tipping elements

# tipping threshold in °C
P_PFTP = np.array([100])  # disable PFTP
P_AMAZ = np.array([100])  # disable AMAZ
P_PFAT = np.array([P_mean['PFAT']])

# Rate - doesn't matter as both elements aren't used here
R_PFTP = np.array([0.1615])
R_AMAZ = np.array([0.1615])

# impact in GtC
K_PFTP = np.array([0])
K_AMAZ = np.array([0])
K_PFAT = np.array([130])

# Feedback in GtC/°C for PFAT
F_PFAT_100 = F_mean['PFAT100']
F_PFAT_300 = F_mean['PFAT300']

# Initial stock
y0_PFTP = 1.32
y0_AMAZ = 0.75
y0_PFAT = 0

# arange parameter sets
TE_params = {'P': np.array([P_PFTP, P_AMAZ, P_PFAT]).T, 'R': np.array([R_PFTP, R_AMAZ]).T,
             'K': np.array([K_PFTP, K_AMAZ, K_PFAT]).T, 'F100': F_PFAT_100, 'F300': F_PFAT_300 - F_PFAT_100,
             'CH4': np.array([CH4_frac['PFTP'], 0, CH4_frac['PFAT']]), 'y0': np.array([y0_PFTP, y0_AMAZ, y0_PFAT])}

# %% get results
result = run_FaIR_intanal(emissions_in=emms,
                          forcing_in=forc,
                          show_run_info=False,
                          TE_params=TE_params)

# %% plot
fig, axes = plt.subplots(3, 1, sharex='col')

PFAT = result['tipping_elements'][ssp]['PFAT']
axes[0].plot(PFAT.C_stock)
axes[1].plot(PFAT.CO2_emm)
axes[2].plot(PFAT.CH4_emm)

plt.show()


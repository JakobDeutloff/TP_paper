"""
Script to calibrate the mean r for PFAT
"""
import numpy as np
from Code.Read_data.read_simplified_RCMIP import *
from Code.FAIR2.fair_runnter_ctem import run_FaIR_ctem
from Code.Read_data.Constants_new import *
from scipy.optimize import minimize
from Code.Calibrate.Functions import *


# %% Define objective function

# Objective Function
def calibrate(rate, Ts_goal):
    TE_params['R'][0][2] = rate
    result = run_FaIR_ctem(emissions_in=emms,
                           forcing_in=forc,
                           show_run_info=False,
                           TE_params=TE_params)

    trg_year = find_trg_year(result['tipping_elements'].xs(('PFAT', 'bool_tip'), level=(1, 2), axis=1))
    thr_year = find_threshold(result['tipping_elements'].xs(('PFAT', 'C_stock'), level=(1, 2), axis=1), 0.995)

    # squared differences between needed Ts and actual Ts
    diff = ((thr_year - trg_year) - Ts_goal) ** 2

    return diff


# %% get emissions

emms = SSP_emms[['ssp585']]
emms.columns = emms.columns.remove_unused_levels()
forc = SSP_forc[['ssp585']]
forc.columns = forc.columns.remove_unused_levels()

# %% get sample of tipping probs, timescales and impacts of tipping elements

# tipping threshold in °C

P_PFTP = np.array([100])  # disable PFTP
P_AMAZ = np.array([100])  # disable AMAZ
P_PFAT = np.array([P_mean['PFAT']])

# Rate - start with initial guess - from old calibration of PFTP
R_PFTP = np.array([0.1615])
R_AMAZ = np.array([0.1615])
R_PFAT = np.array([0.1615])

# impact in GtC
K_PFTP = np.array([0])
K_AMAZ = np.array([0])
K_PFAT = np.array([(K_max['PFAT'] + K_min['PFAT'])/2])

# Feedback of PFAT - not used here because we calibrate it with mean impact
F_PFAT_100 = np.array([0])
F_PFAT_300 = np.array([0])
cal_flag = True

# Initial stock
y0_PFTP = np.array([0])
y0_AMAZ = np.array([0])
y0_PFAT = K_PFAT * 0.005  # defined as 0.5% of impact

# arange parameter sets
TE_params = {'P': np.array([P_PFTP, P_AMAZ, P_PFAT]).T, 'R': np.array([R_PFTP, R_AMAZ, R_PFAT]).T,
             'K': np.array([K_PFTP, K_AMAZ, K_PFAT]).T, 'CH4': np.array([CH4_frac['PFTP'], 0, CH4_frac['PFAT']]),
             'y0': np.array([y0_PFTP, y0_AMAZ, y0_PFAT]), 'F100': F_PFAT_100, 'F300': F_PFAT_300, 'cal_flag': cal_flag}

# %% Fit mean parameter - only Ts importan for PFAT as impact contains timescales implicitly

par_mean = minimize(calibrate, x0=np.array([8.5 / Ts_min['PFAT']]), method='Nelder-Mead',
                    bounds=[[0.01, 0.5]],
                    tol=1e-6, args=[Ts_min['PFAT']])

# %% save mean timescale calibration

rate_ts_PFAT = pd.DataFrame(data=np.array([[Ts_mean['PFAT']], [par_mean.x[0]], [par_mean.fun]]).T,
                            columns=['timescale', 'rate', 'error'], index=[0])
rate_ts_PFAT.to_csv('Data/Calibration/R_Ts_E_PFAT.csv')

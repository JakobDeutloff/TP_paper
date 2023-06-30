"""
Script to calibrate r for AMAZ
"""
import numpy as np
from Code.Read_data.read_simplified_RCMIP import *
from Code.FAIR2.fair.fair_runnter_ctem import run_FaIR_intanal
import matplotlib.pyplot as plt
from Code.Read_data.Constants_new import *
from scipy.optimize import minimize
from Code.Calibrate.Functions import *
import pickle


# %% Define objective function

# Objective Function
def calibrate(rate, Ts_goal):
    TE_params['R'][0][1] = rate
    result = run_FaIR_intanal(emissions_in=emms,
                              forcing_in=forc,
                              show_run_info=False,
                              TE_params=TE_params)

    trg_year = find_trg_year(result['tipping_elements'].xs(('AMAZ', 'bool_tip'), level=(1, 2), axis=1))
    thr_year = find_threshold(result['tipping_elements'].xs(('AMAZ', 'C_stock'), level=(1, 2), axis=1), 0.995)

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
P_PFTP = np.array([100])
P_AMAZ = np.array([P_mean['AMAZ']])  # disable AMAZ
P_PFAT = np.array([100])  # disable PFAT

# Rate - start with initial guess - from old calibration
R_PFTP = np.array([0.1615])
R_AMAZ = np.array([0.1])
R_PFAT = np.array([0.1615])

# impact in GtC
K_PFTP = np.array([0])
K_AMAZ = np.array([(K_max['AMAZ'] + K_min['AMAZ'])/2])
K_PFAT = np.array([0])

# Initial stock
y0_PFTP = np.array([0])
y0_AMAZ = K_AMAZ * 0.005  # defined as 0.5% of impact
y0_PFAT = np.array([0])

# Feedback of PFAT - not used here because we calibrate it with mean impact
F_PFAT_100 = np.array([0])
F_PFAT_300 = np.array([0])
cal_flag = False  # only used for calibration of PFAT


# arange parameter sets
TE_params = {'P': np.array([P_PFTP, P_AMAZ, P_PFAT]).T, 'R': np.array([R_PFTP, R_AMAZ, R_PFAT]).T,
             'K': np.array([K_PFTP, K_AMAZ, K_PFAT]).T,
             'CH4': np.array([CH4_frac['PFTP'], 0, CH4_frac['PFAT']]), 'y0': np.array([y0_PFTP, y0_AMAZ, y0_PFAT]),
             'F100': F_PFAT_100, 'F300': F_PFAT_300, 'cal_flag': cal_flag}

# %% Fit min, mean and max parameters

par_mean = minimize(calibrate, x0=np.array([0.07]), method='Nelder-Mead',
                    bounds=[[0.01, 1]],
                    tol=1e-1, args=[Ts_mean['AMAZ']])

par_min = minimize(calibrate, x0=np.array([0.07]), method='Nelder-Mead',
                   bounds=[[0.01, 1]],
                   tol=1e-6, args=[Ts_min['AMAZ']])

par_max = minimize(calibrate, x0=np.array([0.07]), method='Nelder-Mead',
                   bounds=[[0.01, 1]],
                   tol=1e-6, args=[Ts_max['AMAZ']])

# %% Plot the results for min, max and mean to check if correct

TE_params['R'][0][1] = par_max.x[0]
result_max = run_FaIR_intanal(emissions_in=emms,
                              forcing_in=forc,
                              show_run_info=False,
                              TE_params=TE_params)

TE_params['R'][0][1] = par_mean.x[0]
result_mean = run_FaIR_intanal(emissions_in=emms,
                               forcing_in=forc,
                               show_run_info=False,
                               TE_params=TE_params)

# save mean result for plot calibration
pickle.dump(result_mean, open('Data/Calibration/mean_AMAZ.pkl', 'wb'))

TE_params['R'][0][1] = par_min.x[0]
result_min = run_FaIR_intanal(emissions_in=emms,
                              forcing_in=forc,
                              show_run_info=False,
                              TE_params=TE_params)

fig, axes = plt.subplots(3, 1, sharex='col')

axes[0].plot(result_min['tipping_elements'].xs(('AMAZ', 'C_stock'), level=(1, 2), axis=1))
time = find_threshold(result_min['tipping_elements'].xs(('AMAZ', 'C_stock'), level=(1, 2), axis=1), 0.995) - \
       find_trg_year(result_min['tipping_elements'].xs(('AMAZ', 'bool_tip'), level=(1, 2), axis=1))
axes[0].set_title('Timescale = ' + str(time) + ' Years')

axes[1].plot(result_mean['tipping_elements'].xs(('AMAZ', 'C_stock'), level=(1, 2), axis=1))
time = find_threshold(result_mean['tipping_elements'].xs(('AMAZ', 'C_stock'), level=(1, 2), axis=1), 0.995) - \
       find_trg_year(result_mean['tipping_elements'].xs(('AMAZ', 'bool_tip'), level=(1, 2), axis=1))
axes[1].set_title('Timescale = ' + str(time) + ' Years')

axes[2].plot(result_max['tipping_elements'].xs(('AMAZ', 'C_stock'), level=(1, 2), axis=1))
time = find_threshold(result_max['tipping_elements'].xs(('AMAZ', 'C_stock'), level=(1, 2), axis=1), 0.995) - \
       find_trg_year(result_max['tipping_elements'].xs(('AMAZ', 'bool_tip'), level=(1, 2), axis=1))
axes[2].set_title('Timescale = ' + str(time) + ' Years')

plt.show()

# %% Calibrate whole range which will be needed later for interpolation

ts_coarse = np.array([Ts_min['AMAZ'], Ts_mean['AMAZ'], Ts_max['AMAZ']])
rate_coarse = np.array([par_min.x[0], par_mean.x[0], par_max.x[0]])

ts_fine = np.arange(Ts_min['AMAZ'], Ts_max['AMAZ']+2, 2)
rate_fine = np.zeros(ts_fine.shape)
error = np.zeros(ts_fine.shape)

# Initial run with interpolating the initial guess from min, mean and max estimate
for i in range(len(ts_fine)):
    par = minimize(calibrate, x0=np.array(np.interp(ts_fine[i], ts_coarse, rate_coarse)), method='Nelder-Mead',
                   bounds=[[rate_coarse[2], rate_coarse[0]]],
                   tol=1e-6, args=[ts_fine[i]])

    rate_fine[i] = par.x[0]
    error[i] = par.fun

# %% Second rund to get rid of the errors
for i in range(len(ts_fine)):
    if error[i] > 0:
        par = minimize(calibrate, x0=np.array([rate_coarse[2]]), method='Nelder-Mead',
                       bounds=[[rate_coarse[2], rate_coarse[0]]],
                       tol=1e-6, args=[ts_fine[i]])

        rate_fine[i] = par.x[0]
        error[i] = par.fun

# %% plot error and rate

fig, axes = plt.subplots(2, 1, sharex='col', figsize=(7, 7))

axes[0].plot(ts_fine, rate_fine)
axes[1].plot(ts_fine, error)
plt.show()

# %% Export rates, timescales and errors

rate_ts_AMAZ = pd.DataFrame(data=np.array([ts_fine, rate_fine, error]).T, columns=['timescale', 'rate', 'error'])
rate_ts_AMAZ.to_csv('Data/Calibration/R_Ts_E_AMAZ.csv')

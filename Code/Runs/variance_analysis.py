"""
Script to generate output for variance analysis between random and LH sampling
"""
import numpy as np
from Code.Read_data.read_simplified_RCMIP import *
from Code.FAIR2.fair.fair_runnter_intanal import run_FaIR_intanal, get_gas_parameter_defaults, \
    get_thermal_parameter_defaults
from Code.Read_data.Constants_new import *
from Code.Calibrate.find_distributions import Ts_perf_dists, P_perf_dists, F_beta_dists, K_beta_dists
from scipy.stats import qmc
import pandas as pd

# %% Setup

N_ens = 10  # Number of ensembles
N_sets = 5000  # number of ensemble members

# define elements
elements = ['PFTP', 'AMAZ', 'PFAT']

# Choose SSP
ssp = 'ssp585'

# members
members = np.arange(N_sets)
members = [str(x) for x in members]

# get emissions and forcing
emms = SSP_emms[[ssp]]
emms.columns = emms.columns.remove_unused_levels().droplevel()
emms = pd.concat([emms] * N_sets, axis=1, keys=members)

forc = SSP_forc[[ssp]]
forc.columns = forc.columns.remove_unused_levels().droplevel()
forc = pd.concat([forc] * N_sets, axis=1, keys=members)

# Get gas and thermal params
gas_params = get_gas_parameter_defaults()
idx = gas_params.index
columns = pd.MultiIndex.from_product([members, list(gas_params.columns.levels[1])])
gas_params_ens = pd.DataFrame(index=idx, columns=columns, data=np.tile(gas_params.values, N_sets))

thermal_params = get_thermal_parameter_defaults()
idx = thermal_params.index
columns = pd.MultiIndex.from_product([members, list(thermal_params.columns.levels[1])])
thermal_params_ens = pd.DataFrame(index=idx, columns=columns, data=np.tile(thermal_params.values, N_sets))

# Initialize Output DF
multiindex = pd.MultiIndex.from_product(
    [np.arange(N_ens), np.arange(N_sets), ['PFTP', 'AMAZ', 'PFAT'], ['C_stock', 'CO2_emm', 'CH4_emm', 'bool_tip']])

results_LH = pd.DataFrame(index=emms.index, columns=multiindex)
results_rand = pd.DataFrame(index=emms.index, columns=multiindex)

# %% Run Ensemble

for j in range(N_ens):

    #  get LH sample
    P_sample = {}
    sampler = qmc.LatinHypercube(d=1, seed=None)  # Latin hypercube sampler to generate random numers between [0, 1)

    # tipping threshold in °C
    for element in elements:
        P_sample[element] = P_perf_dists[element].ppf(sampler.random(N_sets)).squeeze()

    # timescale in years
    ts_PFTP = Ts_perf_dists['PFTP'].ppf(sampler.random(N_sets)).squeeze()
    ts_AMAZ = Ts_perf_dists['AMAZ'].ppf(sampler.random(N_sets)).squeeze()

    # Rate - we need this rather than timescales to parameterise the model
    R_PFTP = np.interp(ts_PFTP, r_PFTP.timescale, r_PFTP.rate).squeeze()
    R_AMAZ = np.interp(ts_AMAZ, r_AMAZ.timescale, r_AMAZ.rate).squeeze()
    R_PFAT = np.ones(R_AMAZ.shape) * r_PFAT.rate[0]

    # impact in GtC
    K_PFTP = K_beta_dists['PFTP'].ppf(sampler.random(N_sets)).squeeze()
    K_AMAZ = K_beta_dists['AMAZ'].ppf(sampler.random(N_sets)).squeeze()
    K_PFAT = K_beta_dists['PFAT'].ppf(sampler.random(N_sets)).squeeze()

    # Feedback PFAT
    sample = sampler.random(N_sets)
    F_PFAT_100 = F_beta_dists['PFAT100'].ppf(sample).squeeze()
    F_PFAT_300 = F_beta_dists['PFAT300'].ppf(sample).squeeze()
    cal_flag = False

    # Initial stock
    y0_PFTP = 0.005 * K_PFTP
    y0_AMAZ = 0.005 * K_AMAZ
    y0_PFAT = 0.005 * K_PFAT

    # arange parameter sets
    TE_params_LH = {'P': np.array([P_sample['PFTP'], P_sample['AMAZ'], P_sample['PFAT']]).T,
                    'R': np.array([R_PFTP, R_AMAZ, R_PFAT]).T,
                    'K': np.array([K_PFTP, K_AMAZ, K_PFAT]).T, 'F100': F_PFAT_100, 'F300': F_PFAT_300,
                    'CH4': np.array([CH4_frac['PFTP'], 0, CH4_frac['PFAT']]),
                    'y0': np.array([y0_PFTP, y0_AMAZ, y0_PFAT]), 'cal_flag': cal_flag}

    #  Get random sample
    P_sample = {}
    sampler = np.random  # Random sampler to generate random numers between [0, 1)

    # tipping threshold in °C
    for element in elements:
        P_sample[element] = P_perf_dists[element].ppf(sampler.random(N_sets)).squeeze()

    # timescale in years
    ts_PFTP = Ts_perf_dists['PFTP'].ppf(sampler.random(N_sets)).squeeze()
    ts_AMAZ = Ts_perf_dists['AMAZ'].ppf(sampler.random(N_sets)).squeeze()

    # Rate - we need this rather than timescales to parameterise the model
    R_PFTP = np.interp(ts_PFTP, r_PFTP.timescale, r_PFTP.rate).squeeze()
    R_AMAZ = np.interp(ts_AMAZ, r_AMAZ.timescale, r_AMAZ.rate).squeeze()
    R_PFAT = np.ones(R_AMAZ.shape) * r_PFAT.rate[0]

    # impact in GtC
    K_PFTP = K_beta_dists['PFTP'].ppf(sampler.random(N_sets)).squeeze()
    K_AMAZ = K_beta_dists['AMAZ'].ppf(sampler.random(N_sets)).squeeze()
    K_PFAT = K_beta_dists['PFAT'].ppf(sampler.random(N_sets)).squeeze()

    # Feedback PFAT
    sample = sampler.random(N_sets)
    F_PFAT_100 = F_beta_dists['PFAT100'].ppf(sample).squeeze()
    F_PFAT_300 = F_beta_dists['PFAT300'].ppf(sample).squeeze()
    cal_flag = False

    # Initial stock
    y0_PFTP = 0.005 * K_PFTP
    y0_AMAZ = 0.005 * K_AMAZ
    y0_PFAT = 0.005 * K_PFAT

    # arange parameter sets
    TE_params_rand = {'P': np.array([P_sample['PFTP'], P_sample['AMAZ'], P_sample['PFAT']]).T,
                      'R': np.array([R_PFTP, R_AMAZ, R_PFAT]).T,
                      'K': np.array([K_PFTP, K_AMAZ, K_PFAT]).T, 'F100': F_PFAT_100, 'F300': F_PFAT_300,
                      'CH4': np.array([CH4_frac['PFTP'], 0, CH4_frac['PFAT']]),
                      'y0': np.array([y0_PFTP, y0_AMAZ, y0_PFAT]), 'cal_flag': cal_flag}

    # Run Ensemble members
    res_LH = run_FaIR_intanal(emissions_in=emms,
                              forcing_in=forc,
                              thermal_parameters=thermal_params_ens,
                              gas_parameters=gas_params_ens,
                              show_run_info=False,
                              TE_params=TE_params_LH)

    res_rand = run_FaIR_intanal(emissions_in=emms,
                                forcing_in=forc,
                                thermal_parameters=thermal_params_ens,
                                gas_parameters=gas_params_ens,
                                show_run_info=False,
                                TE_params=TE_params_rand)

    res_LH['tipping_elements'].to_csv('Data/SSP_output/585_var_analysis/' + str(j) + 'LH.csv')
    res_rand['tipping_elements'].to_csv('Data/SSP_output/585_var_analysis/' + str(j) +'Rand.csv')

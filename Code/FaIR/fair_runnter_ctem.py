# An adapted version of FaIR deveolped by Jakob Deutloff which allows for the coupling to CTEM
# as described in Deutloff et al. (2023). The code is largely based on the version of FaIR developed by
# Nick Leach & Stuart Jenkins (2021). The original version of FaIR can be found under
# https://github.com/njleach/FAIR/tree/1945d44c7bcf237307264f8a687f65a70ed0e34f
# distributed under CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)

## import required dependencies ##
import numpy as np
import pandas as pd
import numexpr as ne
import scipy as sp
import os
from tqdm import tqdm
from scipy import integrate


## begin auxiliary functions ##

def return_gas_namelist():
    gas_namelist = ['bc',
                    'bc|aci',
                    'bc|bc_on_snow',
                    'c2f6',
                    'c3f8',
                    'c4f10',
                    'c5f12',
                    'c6f14',
                    'c7f16',
                    'c8f18',
                    'c_c4f8',
                    'carbon_dioxide',
                    'carbon_tetrachloride',
                    'carbon_tetrachloride|o3',
                    'cf4',
                    'cfc11',
                    'cfc113',
                    'cfc113|o3',
                    'cfc114',
                    'cfc114|o3',
                    'cfc115',
                    'cfc115|o3',
                    'cfc11|o3',
                    'cfc12',
                    'cfc12|o3',
                    'ch2cl2',
                    'ch2cl2|o3',
                    'ch3ccl3',
                    'ch3ccl3|o3',
                    'chcl3',
                    'chcl3|o3',
                    'co',
                    'co|o3',
                    'halon1202',
                    'halon1202|o3',
                    'halon1211',
                    'halon1211|o3',
                    'halon1301',
                    'halon1301|o3',
                    'halon2402',
                    'halon2402|o3',
                    'hcfc141b',
                    'hcfc141b|o3',
                    'hcfc142b',
                    'hcfc142b|o3',
                    'hcfc22',
                    'hcfc22|o3',
                    'hfc125',
                    'hfc134a',
                    'hfc143a',
                    'hfc152a',
                    'hfc227ea',
                    'hfc23',
                    'hfc236fa',
                    'hfc245fa',
                    'hfc32',
                    'hfc365mfc',
                    'hfc4310mee',
                    'methane',
                    'methane|strat_h2o',
                    'methane|o3',
                    'methyl_bromide',
                    'methyl_bromide|o3',
                    'methyl_chloride',
                    'methyl_chloride|o3',
                    'nf3',
                    'nh3',
                    'nitrous_oxide',
                    'nitrous_oxide|o3',
                    'nmvoc',
                    'nmvoc|o3',
                    'nox',
                    'nox_avi',
                    'nox_avi|contrails',
                    'nox|o3',
                    'oc',
                    'oc|aci',
                    'sf6',
                    'so2',
                    'so2f2',
                    'so2|aci'
                    ]

    return gas_namelist


def return_empty_emissions(df_to_copy=False, start_year=1765, end_year=2500, timestep=1, scen_names=[0],
                           gases_in=return_gas_namelist()):
    """
    This function returns a dataframe of zeros in the correct format for use in FaIRv2.0.0-alpha. Pass an existing emission/ concentration array to return a corresponding forcing array.
    """

    if type(df_to_copy) == pd.core.frame.DataFrame:
        df = pd.DataFrame(index=df_to_copy.index,
                          columns=pd.MultiIndex.from_product([df_to_copy.columns.levels[0], gases_in],
                                                             names=['Scenario', 'Gas']
                                                             )
                          ).fillna(0).apply(pd.to_numeric)

    else:

        df = pd.DataFrame(index=np.arange(start_year, end_year + 1, timestep),
                          columns=pd.MultiIndex.from_product([scen_names, gases_in],
                                                             names=['Scenario', 'Gas']
                                                             )
                          ).fillna(0).apply(pd.to_numeric)
        df.index.rename('Year', inplace=True)

    return df


def return_empty_forcing(df_to_copy=False, start_year=1765, end_year=2500, timestep=1, scen_names=[0]):
    """
    This function returns a dataframe of zeros in the correct format for use in FaIRv2.0.0-alpha. Pass an existing emission/ concentration array to return a corresponding forcing array.
    """

    if type(df_to_copy) == pd.core.frame.DataFrame:
        df = pd.DataFrame(index=df_to_copy.index,
                          columns=pd.MultiIndex.from_product([df_to_copy.columns.levels[0], ['forcing']],
                                                             names=['Scenario', 'Variable']
                                                             )
                          ).fillna(0).apply(pd.to_numeric)

    else:

        df = pd.DataFrame(index=np.arange(start_year, end_year + 1, timestep),
                          columns=pd.MultiIndex.from_product([scen_names, ['forcing']],
                                                             names=['Scenario', 'Gas']
                                                             )
                          ).fillna(0).apply(pd.to_numeric)
        df.index.rename('Year', inplace=True)

    return df


def input_to_numpy(input_df):
    """
    converts the dataframe input into a numpy array for calculation, dimension order = [name, gas, time/parameter]
    """

    return input_df.values.T.reshape(input_df.columns.levels[0].size, input_df.columns.levels[1].size,
                                     input_df.index.size)


def get_gas_parameter_defaults(choose_gases=return_gas_namelist()):
    """
    This function returns the FaIRv2.0.0-alpha default parameter set for a gas set of your choice. Available gases can be viewed using the return_gas_namelist() function.
    """

    CHOOSE_params = pd.read_csv(
        os.path.join(os.path.dirname(__file__),
                     "util/parameter-sets/Complete_gas_cycle_params.csv"), header=[0, 1],
        index_col=0).reindex(choose_gases, axis=1, level=1)

    return CHOOSE_params


def get_thermal_parameter_defaults(TCR=1.79, RWF=0.552, F_2x=3.759):
    """
    This function returns the FaIRv2.0.0-alpha default climate parameters. Use the kwargs to specify pre-determined climate sensitivities. In both cases, the response timescales d1-3 (and the shortest-timescale coefficient, q1) are set to the central estimate of a CMIP6 inferred distribution constrained with observational warming. The constraint does not significantly affect the central estimates of the prior (ie. raw CMIP6 inference) distribution.
    """

    d1 = 0.903
    d2 = 7.92
    d3 = 355
    q1 = 0.180
    ECS = TCR / RWF

    v1 = (1 - (d1 / 69.66) * (1 - np.exp(-69.66 / d1)))
    v2 = (1 - (d2 / 69.66) * (1 - np.exp(-69.66 / d2)))
    v3 = (1 - (d3 / 69.66) * (1 - np.exp(-69.66 / d3)))

    q3 = (((TCR / F_2x) - q1 * (v1 - v2) - (ECS / F_2x) * v2) / (v3 - v2))
    q2 = (ECS / F_2x - q1 - q3)

    df = pd.DataFrame([[d1, d2, d3], [q1, q2, q3]],
                      index=['d', 'q'],
                      columns=pd.MultiIndex.from_product([['default'], [1, 2, 3]])
                      )

    return df.apply(pd.to_numeric)


## end auxiliary functions ##

## begin underlying model functions ##

def calculate_alpha(G, G_A, T, r, g0, g1, iirf100_max=False):
    iirf100_val = ne.evaluate("abs(r0 + rU * (G-G_A) + rT * T + rA * G_A)",
                              {'r0': r[..., 0], 'rU': r[..., 1], 'rT': r[..., 2], 'rA': r[..., 3], 'G': G, 'G_A': G_A,
                               'T': T})
    if iirf100_max:
        iirf100_val = ne.evaluate("where(iirf100_val>iirf100_max,iirf100_max,iirf100_val)")
    alpha_val = ne.evaluate("g0 * exp(iirf100_val / g1)")

    return alpha_val


def calculate_alpha_const(alpha):
    return alpha[..., 0]


def step_concentration(R_old, G_A_old, E, alpha, a, tau, PI_conc, emis2conc, dt=1):
    decay_rate = ne.evaluate("dt/(alpha*tau)")
    decay_factor = ne.evaluate("exp(-decay_rate)")
    R_new = ne.evaluate(
        "E * a / decay_rate * ( 1. - decay_factor ) + R_old * decay_factor")  # there shouldn't be a dt in the first decay rate
    G_A = ne.evaluate("sum(R_new,axis=4)")
    C = ne.evaluate("PI_conc + emis2conc * (G_A + G_A_old) / 2")

    return C, R_new, G_A


def unstep_concentration(R_old, G_A, alpha, a, tau, PI_conc, emis2conc, dt=1):
    decay_rate = dt / (alpha * tau)
    decay_factor = np.exp(-decay_rate)
    E = ((G_A - np.sum(R_old * decay_factor, axis=-1)) / np.sum(a / decay_rate * (1. - decay_factor), axis=-1))
    R_new = E[..., None] * a * 1 / decay_rate * (1. - decay_factor) + R_old * decay_factor

    return E, R_new


def step_forcing(C, PI_conc, f):
    # if the logarithmic/sqrt term is undefined (ie. C is zero or negative), this contributes zero to the overall
    # forcing. An exception will appear, however.

    logforc = ne.evaluate("f1 * where( (C/PI_conc) <= 0, 0, log(C/PI_conc) )",
                          {'f1': f[..., 0], 'C': C, 'PI_conc': PI_conc})
    linforc = ne.evaluate("f2 * (C - PI_conc)", {'f2': f[..., 1], 'C': C, 'PI_conc': PI_conc})
    sqrtforc = ne.evaluate("f3 * ( (sqrt( where(C<0 ,0 ,C ) ) - sqrt(PI_conc)) )",
                           {'f3': f[..., 2], 'C': C, 'PI_conc': PI_conc})

    RF = logforc + linforc + sqrtforc

    return RF


def step_temperature(S_old, F, q, d, dt=1):
    decay_factor = ne.evaluate("exp(-dt/d)")  # Analytically solved with integrating factor
    S_new = ne.evaluate("q * F * (1 - decay_factor) + S_old * decay_factor")
    T = ne.evaluate("sum( (S_old + S_new)/2, axis=3 )")  # They take linear interpolated value at delta t/2

    return S_new, T


# Tipping Elements model (Bifurcation type) introduced in Deutloff et al. (2023)
def TE_model(S, P, K, R, T, dt):
    # Simplify equation
    a = ne.evaluate("R * (T/P)")
    b = ne.evaluate("(R * T) / (P * K)")
    # Calculate new stock
    S_new = ne.evaluate("((b / a) * (1 - exp(-a * dt)) + exp(-a * dt) * S ** (-1)) ** (-1)")
    # Omit noegative changes - ireversibility and omitting change prior to tipping with K=0
    S_new[S_new < S] = S[S_new < S]
    return S_new


## end underlying model functions ##

## begin model function ##


def run_FaIR_ctem(emissions_in=False,
                  forcing_in=False,
                  gas_parameters=get_gas_parameter_defaults(),
                  thermal_parameters=get_thermal_parameter_defaults(),
                  show_run_info=True,
                  const_alpha=False,
                  start_tip_model=1910,
                  TE_params=False
                  ):
    """
    Runs the development version of the FaIRv2.0 model.

    Model description paper: https://doi.org/10.5194/gmd-2019-390

    Parameters:

    emissions_in (pandas.core.frame.DataFrame strictly with column index as pandas.core.indexes.multi.MultiIndex):
    A pandas DataFrame containing emission data for the desired GHG and aerosol species. The columns most be a
    MultiIndex with [scenarios , species] as the levels. The species must be consistent between scenarios.

    forcing_in (pandas.core.frame.DataFrame strictly with column index as pandas.core.indexes.multi.MultiIndex):
    A pandas DataFrame containing data for aggregated external forcing. The columns most be a MultiIndex with
    [scenarios , forcing] as the levels. Note that the length of the inner column level dimension must be one
    (ie. forcings must be aggregated).

    gas_parameters (pandas.core.frame.DataFrame strictly with column index as pandas.core.indexes.multi.MultiIndex):
    A pandas DataFrame containing the gas cycle parameters for the desired GHG and aerosol species. The columns most be
     a MultiIndex with [parameter set , species] as the levels. The species must be consistent between parameter sets.
     'Indirect' forcings can be specified by adding species with the syntax 'x|y': this means the gas cycle of species
     'x' is used to compute an additional forcing based on the f parameters specified. 'y' designates the name of the
     indirect forcing, such as 'methane|strat_h2o'.

    thermal_parameters (pandas.core.frame.DataFrame strictly with column index as pandas.core.indexes.multi.MultiIndex):
    A pandas DataFrame containing the response parameters used for each box. The columns most be a MultiIndex with
    [parameter set , response box] as the levels. Any number of boxes can be specified by varying the number of
    timescales 'd' and coefficients 'q' supplied.

    show_run_info (bool):
    Specify whether to show information about the current run. Suggest setting to True for normal use, but False if
    optimising parameters or running recursively.

    const_alphy (Bool): Determines whether alpha should be kept constant or is allowed to vary.

    y0 (dict): Initial stocks of tipping points model

    start_tip_model (int): Year in which tipping element model starts running. Negative values of temperature change,
    caused by negative emission values which occur in the historical data should be omitted.
    """

    # Determine the number of scenario runs , parameter sets , gases , integration period, timesteps
    time_index = emissions_in.index
    [(dim_scenario, scen_names), (dim_gas_param, gas_set_names), (dim_thermal_param, thermal_set_names)] = [
        (x.size, list(x)) for x in
        [emissions_in.columns.levels[0], gas_parameters.columns.levels[0], thermal_parameters.columns.levels[0]]]
    gas_names = [x for x in gas_parameters.columns.levels[1] if '|' not in x]
    n_gas = len(gas_names)
    n_forc, forc_names = gas_parameters.columns.levels[1].size, list(gas_parameters.columns.levels[1])
    n_year = time_index.size

    # map the concentrations onto the forcings (ie. so the correct indirect forcing parameters read the correct
    # concentration arrays)
    gas_forc_map = [gas_names.index(forc_names[x].split('|')[0]) for x in np.arange(len(forc_names))]

    names_list = [scen_names, gas_set_names, thermal_set_names, gas_names]
    names_titles = ['Scenario', 'Gas cycle set', 'Thermal set', 'Gas name']
    forc_names_list = [scen_names, gas_set_names, thermal_set_names, forc_names]
    forc_names_titles = ['Scenario', 'Gas cycle set', 'Thermal set', 'Forcing component']

    timestep = np.append(np.diff(time_index), np.diff(time_index)[-1])

    # check if no dimensions are degenerate
    if (set(scen_names) != set(gas_set_names)) & (set(scen_names) != set(thermal_set_names)) & (
            set(gas_set_names) != set(thermal_set_names)):
        gas_shape, gas_slice = [1, dim_gas_param, 1], gas_set_names
        therm_shape, therm_slice = [1, 1, dim_thermal_param], thermal_set_names
    # check if all degenerate
    elif (set(scen_names) == set(gas_set_names)) & (set(scen_names) == set(thermal_set_names)):
        gas_shape, gas_slice = [dim_scenario, 1, 1], scen_names
        therm_shape, therm_slice = [dim_scenario, 1, 1], scen_names
        dim_gas_param = 1
        dim_thermal_param = 1
        [x.pop(1) for x in [names_list, names_titles, forc_names_list, forc_names_titles]]
        [x.pop(1) for x in [names_list, names_titles, forc_names_list, forc_names_titles]]
    # check other possibilities
    else:
        if set(scen_names) == set(gas_set_names):
            gas_shape, gas_slice = [dim_scenario, 1, 1], scen_names
            therm_shape, therm_slice = [1, 1, dim_thermal_param], thermal_set_names
            dim_gas_param = 1
            [x.pop(1) for x in [names_list, names_titles, forc_names_list, forc_names_titles]]
        elif set(scen_names) == set(thermal_set_names):
            gas_shape, gas_slice = [1, dim_gas_param, 1], gas_set_names
            therm_shape, therm_slice = [dim_scenario, 1, 1], scen_names
            dim_thermal_param = 1
            [x.pop(2) for x in [names_list, names_titles, forc_names_list, forc_names_titles]]
        else:
            gas_shape, gas_slice = [1, dim_gas_param, 1], gas_set_names
            therm_shape, therm_slice = [1, dim_gas_param, 1], gas_set_names
            dim_thermal_param = 1
            [x.pop(2) for x in [names_list, names_titles, forc_names_list, forc_names_titles]]

    ## Reindex to align columns:
    emissions = emissions_in.reindex(scen_names, axis=1, level=0).reindex(gas_names, axis=1, level=1).values.T.reshape(
        dim_scenario, 1, 1, n_gas, n_year)

    if forcing_in is False:
        ext_forcing = np.zeros((dim_scenario, 1, 1, 1, n_year))
    else:
        forcing_in = forcing_in.reindex(scen_names, axis=1, level=0)
        ext_forcing = forcing_in.loc[:, (scen_names, slice(None))].values.T.reshape(dim_scenario, 1, 1, 1, n_year)

    gas_cycle_parameters = gas_parameters.reindex(gas_slice, axis=1, level=0).reindex(gas_names, axis=1, level=1)
    thermal_parameters = thermal_parameters.reindex(therm_slice, axis=1, level=0)

    ## get parameter arrays
    a, tau, r, PI_conc, emis2conc = [gas_cycle_parameters.loc[x].values.T.reshape(gas_shape + [n_gas, -1]) for x in
                                     [['a1', 'a2', 'a3', 'a4'], ['tau1', 'tau2', 'tau3', 'tau4'],
                                      ['r0', 'rC', 'rT', 'rA'], 'PI_conc', 'emis2conc']]
    f = gas_parameters.reindex(gas_slice, axis=1, level=0).reindex(forc_names, axis=1, level=1).loc[
        'f1':'f3'].values.T.reshape(gas_shape + [n_forc, -1])
    d, q = [thermal_parameters.loc[x].values.T.reshape(therm_shape + [-1]) for x in ['d', 'q']]

    if show_run_info:
        print('Integrating ' + str(dim_scenario) + ' scenarios, ' + str(
            dim_gas_param) + ' gas cycle parameter sets, ' + str(
            dim_thermal_param) + ' thermal response parameter sets, over ' + str(
            forc_names) + ' forcing agents, between ' + str(time_index[0]) + ' and ' + str(time_index[-1]) + '...',
              flush=True)

    # Dimensions : [scenario, gas params, thermal params, gas, time, (gas/thermal pools)]
    g1 = np.sum(a * tau * (1. - (1. + 100 / tau) * np.exp(-100 / tau)), axis=-1)
    g0 = np.exp(-1 * np.sum(a * tau * (1. - np.exp(-100 / tau)), axis=-1) / g1)

    # Create appropriate shape variable arrays / calculate RF if concentration driven
    C = np.zeros((dim_scenario, dim_gas_param, dim_thermal_param, n_gas, n_year))
    RF = np.zeros((dim_scenario, dim_gas_param, dim_thermal_param, n_forc, n_year))
    T = np.zeros((dim_scenario, dim_gas_param, dim_thermal_param, n_year))
    alpha = np.zeros((dim_scenario, dim_gas_param, dim_thermal_param, n_gas, n_year))
    alpha[..., 0] = calculate_alpha(G=0, G_A=0, T=0, r=r, g0=g0, g1=g1)

    G = np.cumsum(emissions, axis=-1)
    C[..., 0], R, G_A = step_concentration(R_old=0, G_A_old=0, alpha=alpha[..., 0, np.newaxis],
                                           E=emissions[..., 0, np.newaxis], a=a, tau=tau, PI_conc=PI_conc[..., 0],
                                           emis2conc=emis2conc[..., 0], dt=timestep[0])
    RF[..., 0] = step_forcing(C=C[..., gas_forc_map, 0], PI_conc=PI_conc[..., gas_forc_map, 0], f=f)
    S, T[..., 0] = step_temperature(S_old=0, F=np.sum(RF[..., 0], axis=-1)[..., np.newaxis] + ext_forcing[..., 0],
                                    q=q, d=d, dt=timestep[0])

    # Prepare tipping elements model - Definition of tipping element parameters without uncertainty estimates

    # Initial Stocks
    if TE_params:
        dim_tipping_elems = 3
        stock_TE = np.ones((dim_scenario, dim_tipping_elems, n_year))
        stock_TE[:, 0, :] = stock_TE[:, 0, :] * np.expand_dims(TE_params['y0'][0], axis=1)
        stock_TE[:, 1, :] = stock_TE[:, 1, :] * np.expand_dims(TE_params['y0'][1], axis=1)
        stock_TE[:, 2, :] = stock_TE[:, 2, :] * np.expand_dims(TE_params['y0'][2], axis=1)
        emm_TE_CO2 = np.zeros((dim_scenario, dim_tipping_elems, n_year))
        emm_TE_CH4 = np.zeros((dim_scenario, dim_tipping_elems, n_year))
        weight_C = 12.011  # weight of carbon in u
        weight_CH4 = 16  # weight of CH4 in u

        # Implement tipping behaviour, K is set to the given impact once P is crossed
        K = np.zeros(TE_params['K'].shape)
        # Bool array saying whether T crossed P already
        bool = np.full((dim_scenario, dim_tipping_elems, n_year), False)

    # Timestepping
    for t in tqdm(np.arange(1, n_year), unit=' timestep'):

        if time_index[t] == 1901:
            # Calculate pre-industrial temperature
            temp_1850_1900 = T[..., 100:150].mean(axis=3).squeeze(axis=2)

        if TE_params:
            # Calculate emissions from tipping points
            if time_index[t] >= start_tip_model:

                # check if thresholds are crossed
                bool_new = TE_params['P'] <= T[..., t - 1].squeeze(axis=2) - temp_1850_1900
                # update memebers with crossed thresholds to true - can't be reversed
                bool[bool_new, t:] = True
                # set K to given impact if tipped for PFTP and AMAZ
                K[:, :2][bool[:, :2, t]] = TE_params['K'][:, :2][bool[:, :2, t]]

                # calculate K of PFAT
                # assign new variables to make calculation sleeker
                F100 = TE_params['F100'][bool[:, 2, t]]
                F300 = TE_params['F300'][bool[:, 2, t]]
                ti = time_index[t]
                # little trick for calibration where T would be one-dimensional and squeeze doesn't work
                if dim_scenario == 1:
                    temp = T[..., t - 1].reshape(1)[bool[:, 2, t]]
                else:
                    temp = T[..., t - 1].squeeze()[bool[:, 2, t]]
                if time_index[t] <= 2100:
                    K[:, 2][bool[:, 2, t]] = F100 * temp
                elif (time_index[t] > 2100) and (time_index[t] <= 2300):
                    # linear interpolation between feedback parameters
                    K[:, 2][bool[:, 2, t]] = ((F100 * (2300 - ti) + F300 * (ti - 2100))/200) * temp
                elif time_index[t] > 2300:
                    K[:, 2][bool[:, 2, t]] = F300 * temp
                # in case of calibration just use given impact
                if TE_params['cal_flag']:
                    K[:, 2][bool[:, 2, t]] = TE_params['K'][:, 2][bool[:, 2, t]]
                # if calculated K exceeds K_max set K to K_max
                bool_K = K[:, 2] > TE_params['K'][:, 2]
                K[:, 2][bool_K] = TE_params['K'][:, 2][bool_K]

                # Calculate stocks
                stock_TE[:, :, t] = TE_model(stock_TE[:, :, t - 1], TE_params['P'], K, TE_params['R'],
                                             T[..., t - 1].squeeze(axis=2), timestep[t])

                # Add emissions
                idx_CO2 = gas_names.index('carbon_dioxide')
                idx_CH4 = gas_names.index('methane')
                emm_TE_CO2[:, :, t] = (stock_TE[:, :, t] - stock_TE[:, :, t - 1]) * (1 - TE_params['CH4'] / 100)
                emm_TE_CH4[:, :, t] = (stock_TE[:, :, t] - stock_TE[:, :, t - 1]) * (TE_params['CH4'] / 100) * \
                                      (weight_CH4 / weight_C) * 1e3  # conversion factor since CH4 is given in Mt CH4
                emissions[:, 0, 0, idx_CO2, t] += np.sum(emm_TE_CO2[:, :, t], axis=1)
                emissions[:, 0, 0, idx_CH4, t] += np.sum(emm_TE_CH4[:, :, t], axis=1)

        # Calculate new FaIR Timestep
        if const_alpha:
            alpha[..., t] = calculate_alpha_const(alpha)
        else:
            alpha[..., t] = calculate_alpha(G=G[..., t - 1], G_A=G_A, T=np.sum(S, axis=-1)[..., np.newaxis], r=r,
                                            g0=g0, g1=g1)

        C[..., t], R, G_A = step_concentration(R_old=R, G_A_old=G_A, alpha=alpha[..., t, np.newaxis],
                                               E=emissions[..., t, np.newaxis], a=a, tau=tau,
                                               PI_conc=PI_conc[..., 0], emis2conc=emis2conc[..., 0], dt=timestep[t])
        RF[..., t] = step_forcing(C=C[..., gas_forc_map, t], PI_conc=PI_conc[..., gas_forc_map, 0], f=f)
        S, T[..., t] = step_temperature(S_old=S,
                                        F=np.sum(RF[..., t], axis=-1)[..., np.newaxis] + ext_forcing[..., t], q=q,
                                        d=d, dt=timestep[t])

    C_out = pd.DataFrame(np.moveaxis(C, -1, 0).reshape(C.shape[-1], -1), index=time_index,
                         columns=pd.MultiIndex.from_product(names_list, names=names_titles))
    E_out = emissions_in.copy()

    ext_forcing = np.zeros(np.sum(RF, axis=-2)[..., np.newaxis, :].shape) + ext_forcing
    RF = np.concatenate((RF, ext_forcing), axis=-2)
    RF = np.concatenate((RF, np.sum(RF, axis=-2)[..., np.newaxis, :]), axis=-2)

    alpha_out = pd.DataFrame(np.moveaxis(alpha, -1, 0).reshape(alpha.shape[-1], -1), index=time_index,
                             columns=pd.MultiIndex.from_product(names_list, names=names_titles))
    RF_out = pd.DataFrame(np.moveaxis(RF, -1, 0).reshape(RF.shape[-1], -1), index=time_index,
                          columns=pd.MultiIndex.from_product(
                              [x + ['External', 'Total'] * (x == forc_names_list[-1]) for x in forc_names_list],
                              names=forc_names_titles))
    temp = np.moveaxis(T, -1, 0).reshape(T.shape[-1], -1)
    temp_pid = temp[100:150, ...].mean(axis=0)
    T_out = pd.DataFrame(temp - temp_pid, index=time_index,
                         columns=pd.MultiIndex.from_product(names_list[:-1], names=names_titles[:-1]))

    if TE_params:
        # prepare output of tipping elements
        multiindex = pd.MultiIndex.from_product(
            [scen_names, ['PFTP', 'AMAZ', 'PFAT'], ['C_stock', 'CO2_emm', 'CH4_emm', 'bool_tip']])
        data = np.empty((len(time_index), dim_scenario * dim_tipping_elems * 4))
        # first column of each scenario: C stock in GtC
        data[:, np.arange(0, data.shape[1] - 11, 12)] = stock_TE[:, 0, :].T  # PFTP
        data[:, np.arange(4, data.shape[1] - 7, 12)] = stock_TE[:, 1, :].T  # AMAZ
        data[:, np.arange(8, data.shape[1] - 3, 12)] = stock_TE[:, 2, :].T  # PFAT
        # second column of each scenario: CO2 emissions in GtC
        data[:, np.arange(1, data.shape[1] - 10, 12)] = emm_TE_CO2[:, 0, :].T  # PFTP
        data[:, np.arange(5, data.shape[1] - 6, 12)] = emm_TE_CO2[:, 1, :].T  # AMAZ
        data[:, np.arange(9, data.shape[1] - 2, 12)] = emm_TE_CO2[:, 2, :].T  # PFAT
        # third column of each scenario: CH4 emissions in MtC
        data[:, np.arange(2, data.shape[1] - 9, 12)] = emm_TE_CH4[:, 0, :].T  # PFTP
        data[:, np.arange(6, data.shape[1] - 5, 12)] = emm_TE_CH4[:, 1, :].T  # AMAZ
        data[:, np.arange(10, data.shape[1] - 1, 12)] = emm_TE_CH4[:, 2, :].T  # PFAT
        # fourth column of each scenario: bool of tipping
        data[:, np.arange(3, data.shape[1] - 8, 12)] = bool[:, 0, :].T  # PFTP
        data[:, np.arange(7, data.shape[1] - 4, 12)] = bool[:, 1, :].T  # AMAZ
        data[:, np.arange(11, data.shape[1], 12)] = bool[:, 2, :].T  # PFAT
        TE_out = pd.DataFrame(data=data, index=time_index, columns=multiindex)
    else:
        TE_out = pd.DataFrame(data=[])

    out_dict = {'C': C_out,
                'RF': RF_out,
                'T': T_out,
                'alpha': alpha_out,
                'Emissions': E_out,
                'gas_parameters': gas_parameters,
                'thermal parameters': thermal_parameters,
                'tipping_elements': TE_out}

    for axis in [x for x in list(out_dict.keys())[:-2] if type(x) == pd.core.frame.DataFrame]:
        out_dict[axis].index = out_dict[axis].index.rename('Year')

    return out_dict

## end model function ##

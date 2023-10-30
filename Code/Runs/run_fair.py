"""
Script to produce coupled and uncoupled model ensembles
"""

import glob
from Code.Calibrate.find_distributions import Ts_perf_dists, P_perf_dists, F_beta_dists, K_beta_dists
from Code.FaIR.fair_runnter_ctem import run_FaIR_ctem
from Code.Read_data.read_simplified_RCMIP import SSP_emms
from Code.Read_data.Constants_new import CH4_frac
from Code.Read_data.Constants_new import r_PFTP, r_AMAZ, r_PFAT
from scipy.stats import qmc
import pandas as pd
import numpy as np
import pickle


# %% define functions

def load_fair_params(N_sets, from_source, ens_name):
    # Load the categories of the radiative forcing agents that are used later to group radiative forcings
    param_categories = pd.read_csv(r'Data/FaIRv2.0.0-alpha_RF_categories.csv',
                                   index_col=0,
                                   skiprows=1, names=['source', 'category'])
    param_categories.loc['Total'] = 'Total'

    # Get new parameters for gas, thermal and extraforcing if needed
    if not from_source:
        # Select the ensemble members of the unconstrained 1e6 member ensemble for the constrained ensemble
        # The observational datasets of current warming that are used
        datasets_to_use = ['HadCRUT5', 'HadCRUT4', 'NOAA', 'GISTEMP', 'CW', 'BERKELEY']
        # Load the 1e6 probabilities of exceeding the observed warming
        # noinspection PyTypeChecker
        FULL_probabilities = pd.concat(
            [pd.read_hdf(r"C:\Users\jakob\3D Objects\util_data\parameter-sets/perturbed-parameters"
                         r"/FULL_selection_probability-" + x + '.h5') for x in datasets_to_use],
            axis=1, keys=datasets_to_use)
        # Randomly pic the ensemble members for the constrained ensemble with higher probability of picking those
        # agreeing with observed warming trend.
        FULL_ensemble_selection = FULL_probabilities.mean(axis=1) > np.random.random(FULL_probabilities.shape[0])

        # here we randomly select which ensemble members will make up our parameter sets
        CONSTRAINED_ensemble_members = FULL_ensemble_selection.index[FULL_ensemble_selection][
            np.random.choice(FULL_ensemble_selection.sum(), N_sets, replace=False)]

        # Get the parameters for the constrained ensemble members
        # member names, grouping 10000 memebers each to read in parameter data which is stored in seperate files,
        # one for 10000 members each
        ALL_mems = [x.split(r"gas_params")[-1].split('.')[0] for x in
                    glob.glob(
                        r'C:\Users\jakob\3D Objects\util_data\parameter-sets/perturbed-parameters/gas_params/*.h5')]

        CONSTRAINED_thermal_set = []
        CONSTRAINED_gas_set = []
        CONSTRAINED_extforc_sfs = []

        # Read the files for all 10000 member groups and add the ones picked for the constrained ensemble to the
        # constrained parameter sets
        for mem_range in ALL_mems:
            gas_params = pd.read_hdf(
                r'C:\Users\jakob\3D Objects\util_data\parameter-sets/perturbed-parameters/gas_params' + mem_range + '.h5')
            thermal_params = pd.read_hdf(
                r'C:\Users\jakob\3D Objects\util_data\parameter-sets/perturbed-parameters/climresp_params/FULL' + mem_range +
                '.h5')
            extforc_sfs = pd.read_hdf(
                r'C:\Users\jakob\3D Objects\util_data\parameter-sets/perturbed-parameters/extforc_sfs' + mem_range + '.h5')

            # Get index of memebrs that belong to the constrained ensemble
            CONSTRAINED_mems = set(gas_params.columns.levels[0]).intersection(CONSTRAINED_ensemble_members)

            # Add the choosen members to the constrained parameter sets
            # d and q for all three thermal boxes
            CONSTRAINED_thermal_set += [thermal_params.reindex(CONSTRAINED_mems, axis=1, level=0)]
            # gas cycle parameters for all four gas pools and timescale adjustment and radiative forcing parameters
            CONSTRAINED_gas_set += [gas_params.reindex(CONSTRAINED_mems, axis=1, level=0)]
            # Uncertainty in effective radiative forcing which is normally included in f1,2,3 (gas_params) is given
            # explicitly for external forcing here
            CONSTRAINED_extforc_sfs += [extforc_sfs.reindex(CONSTRAINED_mems, axis=1)]

        CONSTRAINED_thermal_set = pd.concat(CONSTRAINED_thermal_set, axis=1)
        CONSTRAINED_gas_set = pd.concat(CONSTRAINED_gas_set, axis=1)
        CONSTRAINED_extforc_sfs = pd.concat(CONSTRAINED_extforc_sfs, axis=1)

        # save newly created parameter sets
        pickle.dump(CONSTRAINED_thermal_set, open('Data/Params/' + ens_name + '/thermal_set.pkl', 'wb'))
        pickle.dump(CONSTRAINED_gas_set, open('Data/Params/' + ens_name + '/gas_set.pkl', 'wb'))
        pickle.dump(CONSTRAINED_extforc_sfs, open('Data/Params/' + ens_name + '/extforc.pkl', 'wb'))

    #  If previously saved parameters should be used, they are loaded
    else:
        ens_nam = 'coupled_ensemble'
        CONSTRAINED_thermal_set = pickle.load(open('Data/Params/' + ens_nam + '/thermal_set.pkl', 'rb'))
        CONSTRAINED_gas_set = pickle.load(open('Data/Params/' + ens_nam + '/gas_set.pkl', 'rb'))
        CONSTRAINED_extforc_sfs = pickle.load(open('Data/Params/' + ens_nam + '/extforc.pkl', 'rb'))

    return CONSTRAINED_thermal_set, CONSTRAINED_gas_set, CONSTRAINED_extforc_sfs, param_categories


def load_TE_params(from_source, elements, N_sets, ens_name):
    # get sample of tipping probs, timescales and impacts of tipping elements

    if not from_source:

        sampler = qmc.LatinHypercube(d=1, seed=None)  # Latin hypercube sampler to generate random numers between [0, 1)

        P_sample = {}
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

        # Feedback PFAT - only correlation included
        sample = sampler.random(N_sets)
        F_PFAT_100 = F_beta_dists['PFAT100'].ppf(sample).squeeze()
        F_PFAT_300 = F_beta_dists['PFAT300'].ppf(sample).squeeze()
        cal_flag = False

        # Initial stock
        y0_PFTP = 0.005 * K_PFTP
        y0_AMAZ = 0.005 * K_AMAZ
        y0_PFAT = 0.005 * K_PFAT

        # arange parameter sets
        TE_params = {'P': np.array([P_sample['PFTP'], P_sample['AMAZ'], P_sample['PFAT']]).T,
                     'R': np.array([R_PFTP, R_AMAZ, R_PFAT]).T,
                     'K': np.array([K_PFTP, K_AMAZ, K_PFAT]).T, 'F100': F_PFAT_100, 'F300': F_PFAT_300,
                     'CH4': np.array([CH4_frac['PFTP'], 0, CH4_frac['PFAT']]),
                     'y0': np.array([y0_PFTP, y0_AMAZ, y0_PFAT]), 'cal_flag': cal_flag}

        pickle.dump(TE_params, open('Data/Params/' + ens_name + '/TE_params.pkl', 'wb'))
        pickle.dump(P_sample, open('Data/Params/' + ens_name + '/P_sample.pkl', 'wb'))

    else:
        ens_nam = 'coupled_ensemble'
        TE_params = pickle.load(open('Data/Params/' + ens_nam + '/TE_params.pkl', 'rb'))
        P_sample = pickle.load(open('Data/Params/' + ens_nam + '/P_sample.pkl', 'rb'))

    return TE_params, P_sample


def run_ssp(ssp, use_TE_model, thermal_set, gas_set, extforc_sfs, TE_params, P_sample, N_sets, param_categories,
            elements, tip_prob_total):
    """
    Function used to run FAiR with respective ssp and calculate output variables from the model results
    :param TE_params:
    :param use_TE_model: bool whether to use TE model
    :param ssp: str of ssp to use
    :return: list of needed variables
    """
    ## load emissions from pre-computed (simplify_RCMIP.py) FAiR emissions
    emms = SSP_emms[[ssp]].loc[1750:]
    emms.columns = emms.columns.remove_unused_levels().droplevel()

    # Load total erf timeseries - not the RMIP dataset used in the default FAiR setting and the calibration of the
    # permafrost module but nearly no difference. Readyly processed forcing used in default conf. of Fair can't be used
    # as we need to distinguish between solar etc. to apply uncertainties
    ssp_erf = pd.read_csv(
        'https://raw.githubusercontent.com/Priestley-Centre/ssp_erf/master/SSPs/ERF_' + ssp + '_1750-2500.csv',
        index_col=0, dtype=float).loc[1750:]

    # Pick land use, volcanic and solar erf timeseries and store the sum of all three for each ensemble member
    # Multiply the erf timeseries with the uncertainty in erf for each member
    extforc_arr = (extforc_sfs.values[None] * ssp_erf[['land_use', 'volcanic', 'solar']].values[..., None])
    extforc_ts = pd.DataFrame(data=extforc_arr.sum(axis=1), index=ssp_erf.index,
                              columns=pd.MultiIndex.from_product([extforc_sfs.columns, ['forcing']]))
    emms = pd.concat([emms] * N_sets, axis=1, keys=extforc_ts.columns.levels[0])

    if use_TE_model:
        result = run_FaIR_ctem(emissions_in=emms,
                               forcing_in=extforc_ts,
                               gas_parameters=gas_set,
                               thermal_parameters=thermal_set,
                               show_run_info=False,
                               TE_params=TE_params)  # set to false to uncouple CTEM

    else:
        result = run_FaIR_ctem(emissions_in=emms,
                               forcing_in=extforc_ts,
                               gas_parameters=gas_set,
                               thermal_parameters=thermal_set,
                               show_run_info=False,
                               TE_params=False)  # set to false to uncouple CTEM

    # Calculate anthropogenic RF
    # RF from Land use cange
    LUC_rf = pd.DataFrame(data=extforc_arr[:, 0, :], index=ssp_erf.index,
                          columns=pd.MultiIndex.from_product([extforc_sfs.columns, ['forcing']])).droplevel(
        axis=1, level=1).stack().sort_index()

    # Define RF as relative to 1850 and sort by focing categories
    rf = (result['RF'] - result['RF'].loc[1850]).stack(level=0).groupby(param_categories.category.to_dict(),
                                                                        axis=1).sum().sort_index()
    # RF from aerosols
    rf['aer'] = rf.aci + rf.ari
    # RF Anthropogenic
    rf['anthro'] = rf.aer + rf.bc_on_snow + rf.carbon_dioxide + rf.contrails + rf.methane + rf.nitrous_oxide + \
                   rf.other_wmghgs + rf.ozone + rf.strat_h2o + LUC_rf

    # Get quantiles
    rf = rf.stack().unstack(level=1).quantile([0.05, 0.166, 0.5, 0.833, 0.95], axis=1).stack().swaplevel(0, 1) \
        .sort_index().T

    T = result['T'].quantile([0.05, 0.166, 0.5, 0.833, 0.95], axis=1).T

    C = result['C'].groupby(level=[1], axis=1).quantile([0.05, 0.166, 0.5, 0.833, 0.95])

    alpha_C = result['alpha'].xs('carbon_dioxide', level='Gas name', axis=1).quantile([0.05, 0.166, 0.5, 0.833, 0.95],
                                                                                      axis=1).T
    alpha_m = result['alpha'].xs('methane', level='Gas name', axis=1).quantile([0.05, 0.166, 0.5, 0.833, 0.95],
                                                                               axis=1).T
    # Emm quantiles
    if use_TE_model:
        Emissions = result['Emissions'].groupby(level=[1], axis=1).quantile([0.05, 0.166, 0.5, 0.833, 0.95])
    else:
        Emissions = result['Emissions'][result['Emissions'].columns[0][0]]

    # TE model quantiles
    if use_TE_model:
        TE_output = result['tipping_elements'].groupby(level=(1, 2), axis=1).quantile([0.05, 0.166, 0.5, 0.833, 0.95])
    else:
        TE_output = pd.DataFrame(data=[])

    # get total tipping probabilities
    members = np.array(list(result['T'].columns.levels[0]))
    for element in elements:
        for i in range(len(members)):
            try:
                tip_year = result['T'].index[np.where(result['T'][members[i]] > P_sample[element][i])[0][0]]
                tip_prob_total.loc[tip_year:, (ssp, element)] += 1
            except:
                pass

    # Extract results from all enseble members for 2050, 2100, 2300, 2500
    years = [2050, 2100, 2300, 2500]
    result_years = {}
    if use_TE_model:
        keys = ['C', 'RF', 'T', 'Emissions', 'tipping_elements']
    else:
        keys = ['C', 'RF', 'T', 'Emissions']

    for key in keys:
        result_years[key] = result[key].loc[years]

    # Extract 50 random members for T and E
    rand_mems = {}
    rand_idx = np.random.randint(0, len(members), 50)
    rand_members = members[rand_idx]
    rand_mems['T'] = result['T'].loc[:, rand_members]
    rand_mems['Emissions'] = result['Emissions'].loc[:, rand_members]
    if use_TE_model:
        rand_mems['tipping_elements'] = result['tipping_elements'].loc[:, rand_members]

    return T, rf, C, alpha_C, alpha_m, Emissions, result['tipping_elements'], TE_output, \
           tip_prob_total, result_years, rand_mems, result['T']


def run_fair_ens(choose_ssps, elements, N_sets, ens_name, use_TE_model, thermal_set, gas_set, extforc_sfs, TE_params,
                 P_sample, param_categories):
    #  Initialize output arrays SSP ensemble
    quant_T = []
    all_T = []
    quant_RF = []
    quant_C = []
    quant_alpha_C = []
    quant_alpha_meth = []
    Emm_ssp = []
    quant_TE = []
    total_TE = []
    res_years = {}
    rand_mems = {}
    multicol = pd.MultiIndex.from_product([choose_ssps, elements])
    tip_prob_total = pd.DataFrame(data=np.zeros((SSP_emms.index.__len__(), choose_ssps.__len__() * elements.__len__())),
                                  index=SSP_emms.index, columns=multicol)

    for ssp in choose_ssps:
        T, RF, C, a_C, a_m, Emm, TE_whole, \
        TE, tip_prob_total, res_50_100_300_500, rand_mem, T_full = run_ssp(ssp, use_TE_model,
                                                                           thermal_set, gas_set,
                                                                           extforc_sfs,
                                                                           TE_params, P_sample, N_sets,
                                                                           param_categories, elements,
                                                                           tip_prob_total)

        quant_T += [T.copy()]
        all_T += [T_full.copy()]
        quant_RF += [RF.copy()]
        quant_C += [C.copy()]
        quant_alpha_C += [a_C.copy()]
        quant_alpha_meth += [a_m.copy()]
        Emm_ssp += [Emm.copy()]
        quant_TE += [TE.copy()]
        total_TE += [TE_whole.copy()]
        res_years[ssp] = res_50_100_300_500
        rand_mems[ssp] = rand_mem

    quant_T = pd.concat(quant_T, axis=1, keys=choose_ssps)
    all_T = pd.concat(all_T, axis=1, keys=choose_ssps)
    quant_RF = pd.concat(quant_RF, axis=1, keys=choose_ssps)
    quant_C = pd.concat(quant_C, axis=1, keys=choose_ssps)
    quant_alpha_C = pd.concat(quant_alpha_C, axis=1, keys=choose_ssps)
    quant_alpha_meth = pd.concat(quant_alpha_meth, axis=1, keys=choose_ssps)
    Emm_ssp = pd.concat(Emm_ssp, axis=1, keys=choose_ssps)
    quant_TE = pd.concat(quant_TE, axis=1, keys=choose_ssps)
    total_TE = pd.concat(total_TE, axis=1, keys=choose_ssps)
    tip_prob_total = tip_prob_total / N_sets

    #  save results
    quant_RF.to_csv('Data/model_output/' + ens_name + '/RF.csv')
    quant_T.to_csv('Data/model_output/' + ens_name + '/T.csv')
    quant_C.to_csv('Data/model_output/' + ens_name + '/C.csv')
    quant_alpha_C.to_csv('Data/model_output/' + ens_name + '/alpha_C.csv')
    quant_alpha_meth.to_csv('Data/model_output/' + ens_name + '/alpha_m.csv')
    tip_prob_total.to_csv('Data/model_output/' + ens_name + '/tip_prob_total.csv')
    Emm_ssp.to_csv('Data/model_output/' + ens_name + '/Emm.csv')
    quant_TE.to_csv('Data/model_output/' + ens_name + '/TE.csv')
    pickle.dump(total_TE, open('Data/model_output/' + ens_name + '/TE_total.pkl', 'wb'))
    pickle.dump(res_years, open('Data/model_output/' + ens_name + '/result_years.pkl', 'wb'))
    pickle.dump(rand_mems, open('Data/model_output/' + ens_name + '/rand_mems.pkl', 'wb'))
    pickle.dump(all_T, open('Data/model_output/' + ens_name + '/full_T.pkl', 'wb'))


# main
def main():
    # Set-Up

    # Set number of enesemble members
    N_sets = 5000

    # Specify if saved parameters should be used
    from_source = True

    # Set Ens name - also used to load and save parameters
    ens_name = 'uncoupled_ensemble'

    # If TE model should be coupled
    use_TE_model = True

    # Choose SSPs
    choose_ssps = ['ssp126', 'ssp245', 'ssp370', 'ssp434', 'ssp460', 'ssp585']

    # Choose elements for probabilities
    elements = ['PFTP', 'AMOC', 'GRIS', 'WAIS', 'AMAZ', 'BORF', 'AWSI', 'EAIS', 'EASB', 'GLCR', 'LABC', 'TUND', 'PFAT',
                'BARI', 'REEF', 'SAHL']

    # Get fair params
    thermal_set, gas_set, extforc_sfs, param_categories = load_fair_params(N_sets, from_source, ens_name)
    # Get TE model params
    TE_params, P_sample = load_TE_params(from_source, elements, N_sets, ens_name)
    # Run fair ensemble
    run_fair_ens(choose_ssps, elements, N_sets, ens_name, use_TE_model, thermal_set, gas_set, extforc_sfs, TE_params,
                 P_sample, param_categories)


if __name__ == '__main__':
    main()

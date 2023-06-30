"""
Script to run nine additional coupled ensembles for the histogramms of T and the probabilities of exceeding a certain
temperature target
"""


from Code.Read_data.read_simplified_RCMIP import SSP_emms
import pandas as pd
import numpy as np
import pickle
from Code.Runs.run_fair import load_fair_params, load_TE_params, run_ssp


def run_fair_ens_T(choose_ssps, elements, N_sets, ens_name, use_TE_model, thermal_set, gas_set, extforc_sfs, TE_params,
                   P_sample, param_categories, N_ens):
    
    #  Only res year is needed. tip_prob_total given as well to avoid changing the functions
    res_years = {}
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

        res_years[ssp] = res_50_100_300_500

    pickle.dump(res_years, open('Data/SSP_output/' + ens_name + '/result_years' + str(N_ens) + '.pkl', 'wb'))


# main
def main():
    # Set-Up

    # Set numer of ensembles
    N_ens = 9

    # Set number of enesemble members
    N_sets = 5000

    # Specify if saved parameters should be used
    from_source = False

    # Set Ens name - also used to load and save parameters
    ens_name = '5000_tip_T'

    # If TE model should be coupled
    use_TE_model = True

    # Choose SSPs
    choose_ssps = ['ssp119', 'ssp245']

    # Choose elements for probabilities
    elements = ['PFTP', 'AMOC', 'GRIS', 'WAIS', 'AMAZ', 'BORF', 'AWSI', 'EAIS', 'EASB', 'GLCR', 'LABC', 'TUND', 'PFAT',
                'BARI', 'REEF', 'SAHL']

    for i in range(N_ens):
        # Get fair params
        thermal_set, gas_set, extforc_sfs, param_categories = load_fair_params(N_sets, from_source, ens_name)
        # Get TE model params
        TE_params, P_sample = load_TE_params(from_source, elements, N_sets, ens_name)
        # Run fair ensemble
        run_fair_ens_T(choose_ssps, elements, N_sets, ens_name, use_TE_model, thermal_set, gas_set, extforc_sfs, TE_params,
                     P_sample, param_categories, i)


if __name__ == '__main__':
    main()

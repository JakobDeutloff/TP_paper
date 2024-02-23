"""
This script reads output from the model ensembles, either coupled or uncoupled
"""
import pandas as pd
import pickle
import numpy as np


def read_SSP_outout(ens_name):
    T = pd.read_csv('Data/model_output/' + ens_name + '/T.csv', index_col=0, parse_dates=[0], header=[0, 1])
    T.columns.names = ['scenario', 'percentile']
    RF = pd.read_csv('Data/model_output/' + ens_name + '/RF.csv', index_col=0, parse_dates=[0], header=[0, 1, 2])
    RF.columns.names = ['scenario', 'forcing component', 'percentile']
    C = pd.read_csv('Data/model_output/' + ens_name + '/C.csv', index_col=0, parse_dates=[0], header=[0, 1, 2])
    C.columns.names = ['scenario', 'gas', 'percentile']
    alpha_C = pd.read_csv('Data/model_output/' + ens_name + '/alpha_C.csv', index_col=0, parse_dates=[0], header=[0, 1])
    alpha_C.columns.names = ['scenario', 'percentile']
    alpha_m = pd.read_csv('Data/model_output/' + ens_name + '/alpha_m.csv', index_col=0, parse_dates=[0], header=[0, 1])
    alpha_m.columns.names = ['scenario', 'percentile']
    Emm = pd.read_csv('Data/model_output/' + ens_name + '/Emm.csv', index_col=0, parse_dates=[0], header=[0, 1])
    Emm.columns.names = ['scenario', 'variable']
    if np.isnan(Emm['ssp126']['bc'].iloc[1]):  # Emm output is different for coupled runs as percentiles are calculated
        Emm = pd.read_csv('Data/model_output/' + ens_name + '/Emm.csv', index_col=0, parse_dates=[0], header=[0, 1, 2])
        Emm.columns.names = ['scenario', 'variable', 'percentile']
    tip_probs_total = pd.read_csv('Data/model_output/' + ens_name + '/tip_prob_total.csv', index_col=0, parse_dates=[0],
                                  header=[0, 1])
    tip_probs_total.columns.names = ['scenario', 'element']

    return {'T': T, 'RF': RF, 'C': C, 'Alpha_C': alpha_C, 'Alpha_m': alpha_m, 'Emm': Emm,
            'Tip_probs_total': tip_probs_total}


def read_probs(ens_name):
    tip_years_const = pd.read_csv('Data/model_output/' + ens_name + '/Y_tip.csv', index_col=[0], header=[0, 1])
    tip_emms_const = pd.read_csv('Data/model_output/' + ens_name + '/Cum_C_tip.csv', index_col=[0], header=[0, 1])
    cum_TEs = pd.read_csv('Data/model_output/' + ens_name + '/Cum_TEs.csv', index_col=[0], header=[0, 1])
    slow_tip_probs = pd.read_csv('Data/model_output/' + ens_name + '/slow_tip_probs.csv', index_col=[0], header=[0, 1])
    return {'Years': tip_years_const, 'Emms': tip_emms_const, 'cumTEs': cum_TEs, 'slow_tip_probs': slow_tip_probs}


def read_TE_output(ens_name):
    TE = pd.read_csv('Data/model_output/' + ens_name + '/TE.csv', index_col=0, parse_dates=[0],
                     header=[0, 1, 2, 3])
    TE_total = pickle.load(open('Data/model_output/' + ens_name + '/TE_total.pkl', 'rb'))
    return {'TE': TE, 'TE_total': TE_total}


def read_res_years(ens_name):
    res_years = pickle.load(open('Data/model_output/' + ens_name + '/result_years.pkl', 'rb'))
    return res_years


def read_rand_mems(ens_name):
    rand_mems = pickle.load(open('Data/model_output/' + ens_name + '/rand_mems.pkl', 'rb'))
    return rand_mems


def read_T_ens():
    ensembles = []
    for j in range(9):
        ens = pickle.load(open('Data/model_output/5000_tip_T/result_years' + str(j) + '.pkl', 'rb'))
        ensembles += [ens]
    return ensembles

def read_full_T(ens_name):
    T = pickle.load(open('Data/model_output/' + ens_name + '/full_T.pkl', 'rb'))
    return T

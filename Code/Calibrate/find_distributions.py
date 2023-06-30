"""
Script for determining and plotting the probability functions of tipping threshold temperatures and timescales
"""
import matplotlib.pyplot as plt
import numpy as np
from Code.Read_data.Constants_new import par_P_lognorm_dist, par_P_perf_dist, par_Ts_perf_dist, K_max, K_min, \
    F_min, F_max
from scipy import stats

elements = ['PFTP', 'AMOC', 'GRIS', 'WAIS', 'AMAZ', 'BORF', 'AWSI', 'EAIS', 'EASB', 'GLCR', 'LABC', 'TUND', 'PFAT',
            'BARI', 'REEF', 'SAHL']
ts_elements = ['PFTP', 'AMAZ']

# %% find distributions of tipping points including only as lognorm distributions
P_log_norm_dist = {}
for element in elements:
    sample = np.exp(stats.norm.rvs(loc=par_P_lognorm_dist[element][0], scale=par_P_lognorm_dist[element][1], size=10000))
    params = stats.lognorm.fit(sample)
    P_log_norm_dist[element] = stats.lognorm(s=params[0], loc=params[1], scale=[params[2]])


# %% find distributions of tipping points including all probability distributions

def gompertz_pdf(x, b=0.84, a=0.01):
    """
    Function used to find scipy.stats.gompertz()
    :param x: x values
    :param b: scale
    :param a: shape
    :return: gompertz of x
    """
    return (a/b) * np.exp(x/b) * np.exp(-a*(np.exp(x/b) - 1))

def gompertz_cdf(x, b=0.84, a=0.01):

    """
    :param x: x-vector
    :param b: scale
    :param a: shape
    :return: gompertz cdf
    """

    return 1 - np.exp(-a * (np.exp(x/b) - 1))


P_perf_dists = {}
for element in elements:

    # Triangular PDFs
    if par_P_perf_dist[element]['name'] == 'triang':
        scale = par_P_perf_dist[element]['max'] - par_P_perf_dist[element]['min']
        shape = (par_P_perf_dist[element]['mode'] - par_P_perf_dist[element]['min']) / scale
        loc = par_P_perf_dist[element]['min']
        P_perf_dists[element] = stats.triang(c=shape, loc=loc, scale=scale)

    # log-normal PDFs
    if par_P_perf_dist[element]['name'] == 'lnorm':
        params = par_P_perf_dist[element]
        P_perf_dists[element] = stats.lognorm(s=params['sdlog'], scale=np.exp(params['meanlog']))

    # Normal PDFs
    if par_P_perf_dist[element]['name'] == 'norm':
        P_perf_dists[element] = stats.norm(loc=par_P_perf_dist[element]['mean'], scale=par_P_perf_dist[element]['sd'])

    # Gompertz PDFs
    if par_P_perf_dist[element]['name'] == 'gompertz':
        xx = np.linspace(0, 10, 1000)
        values = gompertz_pdf(xx, par_P_perf_dist[element]['scale'], par_P_perf_dist[element]['shape'])
        sample_gom = np.empty([0])
        for i in range(len(xx)):
            sample_gom = np.append(sample_gom, [xx[i]] * int(np.round(values[i] * 1e3)))

        params = stats.gompertz.fit(sample_gom)
        P_perf_dists[element] = stats.gompertz(c=params[0], loc=params[1], scale=params[2])

#%% Find timescale distributions

Ts_perf_dists = {}
for element in ts_elements:

    # log-normal PDFs
    if par_Ts_perf_dist[element]['name'] == 'lnorm':
        pars = par_Ts_perf_dist[element]
        Ts_perf_dists[element] = stats.lognorm(s=pars['sdlog'], scale=np.exp(pars['meanlog']))

# %% Find impact distributions
K_beta_dists = {}

for element in ['PFTP', 'AMAZ', 'PFAT']:
    sample = np.arange(K_min[element], K_max[element]+1)
    params = stats.beta.fit(sample)
    K_beta_dists[element] = stats.beta(a=1, b=1, loc=params[2], scale=params[3])

# %% Find Feedback distributions
F_beta_dists = {}

for element in ['PFAT100', 'PFAT300']:
    sample = np.arange(F_min[element], F_max[element] + 1)
    params = stats.beta.fit(sample)
    F_beta_dists[element] = stats.beta(a=1, b=1, loc=params[2], scale=params[3])





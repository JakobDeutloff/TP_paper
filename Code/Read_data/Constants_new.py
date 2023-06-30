"""
Definitions of constants used within this work
"""

import pandas as pd

# Impact in GtC - for PFAT and PFGT maximum impact
K_min = {'AMAZ': 30, 'PFTP': 125, 'PFAT': 65}
K_max = {'AMAZ': 75, 'PFTP': 250, 'PFAT': 195}
# conversion factor of 1.4 could be used to translate GtC to GtCe for PFAT not done yet, methane emissions could also
# be included explicitly

# Feedback after tipping  in GtC/°C for year 2100 and 2300
# Central estimates from Daves paper
F_mean = {'PFAT100': 10, 'PFGT100': 20, 'PFAT300': 25, 'PFGT300': 50}

# Range from Daves paper abrupt thaw adds to gradual thaw
range_PFAT = [0.25, 0.75]

# max and min feedbacks for PFAT
F_min = {'PFAT100': range_PFAT[0]*F_mean['PFGT100'], 'PFAT300': range_PFAT[0]*F_mean['PFGT300']}
F_max = {'PFAT100': range_PFAT[1]*F_mean['PFGT100'], 'PFAT300': range_PFAT[1]*F_mean['PFGT300']}

# CH4 fraction in %
CH4_frac = {'PFTP': 2.3, 'PFAT': 20, 'PFGT': 2.3}

# Tipping Timescales
Ts_min = {'AMAZ': 50, 'PFTP': 10, 'PFAT': 100, 'PFGT': 100}
Ts_mean = {'AMAZ': 100, 'PFTP': 50, 'PFAT': 200, 'PFGT': 200}
Ts_max = {'AMAZ': 200, 'PFTP': 300, 'PFAT': 300, 'PFGT': 300}

# Rate
r_PFTP = pd.read_csv('Data/Calibration/R_Ts_E_PFTP.csv', header=[0], index_col=[0])
r_AMAZ = pd.read_csv('Data/Calibration/R_Ts_E_AMAZ.csv', header=[0], index_col=[0])
r_PFAT = pd.read_csv('Data/Calibration/R_Ts_E_PFAT.csv', header=[0], index_col=[0])

# mean tip threshold values
P_mean = {'PFTP': 4.0, 'AMOC': 4.0, 'GRIS': 1.5, 'WAIS': 1.5, 'AMAZ': 3.5, 'BORF': 4.0, 'AWSI': 6.3, 'EAIS': 7.5,
          'EASB': 3.0, 'GLCR': 2.0, 'LABC': 1.8, 'TUND': 4.0, 'PFAT': 1.5, 'BARI': 1.6, 'REEF': 1.5, 'SAHL': 2.8}

# max tip threshold
P_max = {'PFTP': 6.0, 'AMOC': 8.0, 'GRIS': 3.0, 'WAIS': 3.0, 'AMAZ': 6.0, 'BORF': 5.0, 'AWSI': 8.7, 'EAIS': 10.0,
         'EASB': 6.0, 'GLCR': 3.0, 'LABC': 3.8, 'TUND': 7.2, 'PFAT': 2.3, 'BARI': 1.7, 'REEF': 2.0, 'SAHL': 3.5}

# min tip threshold
P_min = {'PFTP': 3.0, 'AMOC': 1.4, 'GRIS': 0.8, 'WAIS': 1.0, 'AMAZ': 2.0, 'BORF': 1.4, 'AWSI': 4.5, 'EAIS': 5.0,
         'EASB': 2.0, 'GLCR': 1.5, 'LABC': 1.1, 'TUND': 1.5, 'PFAT': 1.0, 'BARI': 1.5, 'REEF': 1.0, 'SAHL': 2.0}

# lognorm parameters from R script
par_P_lognorm_dist = {'PFTP': [1.3905131, 0.1936294], 'AMOC': [1.3775830, 0.4620041], 'GRIS': [0.4091626, 0.3995689],
                      'WAIS': [0.4077933, 0.2591395], 'AMAZ': [1.251553, 0.333685], 'BORF': [1.386294, 0.135659],
                      'AWSI': [1.839747, 0.200198], 'EAIS': [2.0106831, 0.1936288], 'EASB': [1.1009393, 0.2591508],
                      'GLCR': [0.6973671, 0.1936288], 'LABC': [0.5939343, 0.3281060], 'TUND': [1.3821622, 0.3790672],
                      'PFAT': [0.4067502, 0.2527818], 'BARI': [0.46977524, 0.03795471], 'REEF': [0.4012453, 0.1936286],
                      'SAHL': [1.0267692, 0.1488886]}

# Parameters for fitted probability distributions of tipping point - all distributions
par_P_perf_dist = {'PFTP': {'name': 'triang', 'min': 2.82, 'mode': 2.97, 'max': 6.76},
                   'AMOC': {'name': 'triang', 'min': 0.36, 'mode': 2.64, 'max': 9.85},
                   'GRIS': {'name': 'lnorm', 'meanlog': 0.41, 'sdlog': 0.40},
                   'WAIS': {'name': 'lnorm', 'meanlog': 0.41, 'sdlog': 0.26},
                   'AMAZ': {'name': 'lnorm', 'meanlog': 1.25, 'sdlog': 0.33},
                   'BORF': {'name': 'gompertz', 'shape': 0.007680602, 'scale': 0.841036917},
                   'AWSI': {'name': 'lnorm', 'meanlog': 1.84, 'sdlog': 0.20},
                   'EAIS': {'name': 'norm', 'mean': 7.50, 'sd': 1.52},
                   'EASB': {'name': 'lnorm', 'meanlog': 1.10, 'sdlog': 0.26},
                   'GLCR': {'name': 'lnorm', 'meanlog': 0.70, 'sdlog': 0.19},
                   'LABC': {'name': 'lnorm', 'meanlog': 0.59, 'sdlog': 0.33},
                   'TUND': {'name': 'triang', 'min': 0.38, 'mode': 3.40, 'max': 8.68},
                   'PFAT': {'name': 'lnorm', 'meanlog': 0.41, 'sdlog': 0.25},
                   'BARI': {'name': 'norm', 'mean': 1.60, 'sd': 0.06},
                   'REEF': {'name': 'norm', 'mean': 1.5, 'sd': 0.3},
                   'SAHL': {'name': 'triang', 'min': 1.63, 'mode': 2.88, 'max': 3.82}}

# Parameters for fitted probability distributions of timescales - all distributions
par_Ts_perf_dist = {'PFTP': {'name': 'lnorm', 'meanlog': 3.92, 'sdlog': 1.03},
                    'AMAZ': {'name': 'lnorm', 'meanlog': 4.6051879, 'sdlog': 0.4214717},
                    'PFAT': {'name': 'logistic', 'location': 200, 'scale': 33.96}}

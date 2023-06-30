import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from Code.Read_data.Read_SSP_output import T_const
from Code.Read_data.Constants_new import R_mean, P_mean, Im, Sin_mean
from scipy import integrate

# %%
def Model_PFTP(t, sPFTP, changeGMST, P, Im, R):

    sGMST = np.interp(t, [0, 1], changeGMST)
    dPFTP = R * (sGMST / P) * sPFTP * (1 - sPFTP / Im)

    return dPFTP

def Analytical_PFTP(sPFTP, P, Im, R, T, dt):
    a = R * (T/P)
    b = (R*T)/(P*Im)

    return ((b/a) * (1 - np.exp(-a*dt)) + np.exp(-a*dt) * sPFTP**(-1))**(-1)

ssp = 'ssp585'

start_year = 1910
end_year = 2500
timesteps = end_year - start_year
dt = 1  # timestep in years

PFTP_num = pd.DataFrame(index=np.arange(start_year, end_year+1), columns=['cummC'])
PFTP_an = pd.DataFrame(index=np.arange(start_year, end_year+1), columns=['cummC'])
y = [Sin_mean['PFTP']]
PFTP_num.loc[start_year] = y[0]
PFTP_an.loc[start_year] = y[0]

for i in np.arange(start_year+1, end_year+1):
    changeGMST = [T_const[ssp]['0.5'][i-1], T_const[ssp]['0.5'][i]]
    res = integrate.solve_ivp(Model_PFTP, t_span=[0, 1], y0=PFTP_num.loc[i-1], t_eval=[1],
                              args=[changeGMST, P_mean['PFTP'], Im['PFTP'], R_mean['PFTP']])
    PFTP_num.loc[i] = res.y.squeeze()
    PFTP_an.loc[i] = Analytical_PFTP(PFTP_an.loc[i-1], P_mean['PFTP'], Im['PFTP'], R_mean['PFTP'],
                                     T_const[ssp]['0.5'][i-1], dt)


# %% Test plot
fig, ax = plt.subplots()
ax.plot(PFTP_num, color='b')
ax.plot(PFTP_an, color='r', linestyle='--')
plt.show()

#%%

F100 = 20  # in GtC/°C
F300 = 30  # in GtC/°C
PFGT = []  # Cumulative Emm in GtC

def PFGT(T, x):
    return F100 * T[x] + F300 * T[x-200]






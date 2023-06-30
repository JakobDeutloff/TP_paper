# scripts to help tuning parameters
# Nick Leach (2021)

import numpy as np
import pandas as pd
from fair import calculate_alpha,unstep_concentration


def invert_carbon_cycle_prescribed_T(C,T,a,tau,r,PI_conc,emis2conc):
    
    g1 = np.sum( a * tau * ( 1. - ( 1. + 100/tau ) * np.exp(-100/tau) ), axis=-1 )
    g0 = np.exp( -1 * np.sum( a * tau * ( 1. - np.exp(-100/tau) ) , axis=-1) / g1 ) 
    
    diagnosed_emissions = np.zeros(C.size)
    alpha = np.zeros(C.size)
    G_A = (np.array([np.mean(C[i:i+2]) for i in np.arange(C.size)])-PI_conc)/emis2conc
    G_A[-1]=2*G_A[-1]-G_A[-2]
    
    alpha[0] = calculate_alpha(G=0,G_A=0,T=0,r=r,g0=g0,g1=g1)
    diagnosed_emissions[0],R = unstep_concentration(R_old=0,G_A=G_A[0],alpha=alpha[0,np.newaxis],a=a,tau=tau,PI_conc=PI_conc,emis2conc=emis2conc)
    for t in np.arange(1,C.size):
        G = np.sum(diagnosed_emissions)
        alpha[t] = calculate_alpha(G=G,G_A=G_A[t-1],T=T[t-1],r=r,g0=g0,g1=g1)
        diagnosed_emissions[t],R = unstep_concentration(R_old=R,G_A=G_A[t],alpha=alpha[t,np.newaxis],a=a,tau=tau,PI_conc=PI_conc,emis2conc=emis2conc)
            
    return pd.DataFrame(index=np.arange(C.size),data=np.vstack((diagnosed_emissions,alpha)).T,columns=['emms','alpha'])

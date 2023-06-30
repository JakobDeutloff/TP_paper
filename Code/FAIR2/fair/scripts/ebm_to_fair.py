# scripts to convert EBM parameters to FaIR compatible parameters
# Nick Leach (2021)

import numpy as np
import scipy as sp
import pandas as pd

## builds the state-space model described by the input parameters
def BuildMat(params):
    """
    Creates matrices describing the state-space EBM formulation used by Cummins et al., 2020. 
    params is an array of EBM parameters in the following order:
    [gamma, C1, C2, C3, kap1, kap2, kap3, epsilon, stds, stdx, F_4x]
    Parameter descriptions are as DOI: 10.1175/JCLI-D-19-0589.1
    """

    A = np.array([[-1*params[0],0,0,0],\
                     [1/params[1],-1*(params[4]+params[5])/params[1],params[5]/params[1],0],\
                     [0,params[5]/params[2],-1*(params[5]+params[7]*params[6])/params[2],params[7]*params[6]/params[2]],\
                     [0,0,params[6]/params[3],-1*params[6]/params[3]]])
    k = A.shape[0]
    b = np.array([params[0],0,0,0]).T
    Q = np.zeros((4,4))
    Q[0,0] = params[8]**2
    Q[1,1] = (params[9]/params[1])**2
    A_d = sp.linalg.expm(A)
    b_d = sp.linalg.solve(A, (A_d - np.identity(k)) @ b)
    ## use Van Loan (1978) to compute the matrix exponential
    H = np.zeros((k*2,k*2))
    H[:k,:k] = -A
    H[:k,k:] = Q
    H[k:,k:] = A.T
    G = sp.linalg.expm(H)
    Q_d = G[k:,k:].T @ G[:k,k:]
    C_d = np.array([[0,1,0,0],\
                   [1,-1*params[4],(1-params[7])*params[6],-1*(1-params[7])*params[6]]])
    
    return A,b,Q,A_d,b_d,Q_d,C_d

## converts the input EBM parameters into IR parameters
def EBM_to_FaIR(params):
    
    """
    params is an array of EBM parameters in the following order:
    [gamma, C1, C2, C3, kap1, kap2, kap3, epsilon, stds, stdx, F_4x]
    """

    A,b,Q,A_d,b_d,Q_d,C_d = BuildMat(params)
    eigval,eigvec = np.linalg.eig(A[1:,1:])
    tau = -1/eigval
    q = tau * ( eigvec[0,:] * np.linalg.inv(eigvec)[:,0] ) / params[1]
    
    return pd.DataFrame([tau,q],index=['d','q'],columns=[1,2,3])
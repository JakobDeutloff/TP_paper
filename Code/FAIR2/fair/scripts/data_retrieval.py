# scripts to retrieve data used in Leach et al., 2020
# Nick Leach (2021)

import pandas as pd
from pathlib import Path

# RCMIP data:

## read in the input maps:

FaIR_to_RCMIP_map = pd.read_csv(Path(__file__).parent / "../util/FaIRv2.0.0-alpha_RCMIP_inputmap.csv",index_col=0)
FaIR_to_RCMIP_map_forc = pd.read_csv(Path(__file__).parent / "../util/FaIRv2.0.0-alpha_RCMIP_inputmap.csv",index_col=0)


## pre-load the data upon import 
def get_RCMIP_data():
    
    RCMIP_concs = pd.read_csv('https://rcmip-protocols-au.s3-ap-southeast-2.amazonaws.com/v5.1.0/rcmip-concentrations-annual-means-v5-1-0.csv').set_index(['Region','Scenario','Variable'])
    
    RCMIP_emms = pd.read_csv('https://rcmip-protocols-au.s3-ap-southeast-2.amazonaws.com/v5.1.0/rcmip-emissions-annual-means-v5-1-0.csv').set_index(['Region','Scenario','Variable'])
    
    RCMIP_forc = pd.read_csv('https://rcmip-protocols-au.s3-ap-southeast-2.amazonaws.com/v5.1.0/rcmip-radiative-forcing-annual-means-v5-1-0.csv').set_index(['Region','Scenario','Variable'])
    
    return RCMIP_concs, RCMIP_emms, RCMIP_forc

RCMIP_concs, RCMIP_emms, RCMIP_forc = get_RCMIP_data()

## functions to generate FaIRv2.0.0-alpha input for given scenarios ##
def RCMIP_to_FaIR_input_emms(scenario):
    _out = RCMIP_emms.loc[('World',scenario,FaIR_to_RCMIP_map.RCMIP_emms_key)].droplevel(level=[0,1]).iloc[:,4:]
    _out_index = _out.index.map({v: k for k, v in FaIR_to_RCMIP_map.RCMIP_emms_key.to_dict().items()})
    _out.set_index(_out_index,inplace=True)
    _out = _out.T.mul(FaIR_to_RCMIP_map.RCMIP_emms_scaling)
    _out.index = _out.index.astype(int)
    return _out.apply(pd.to_numeric)

def get_RCMIP_forc(scenario,drivers=['Effective Radiative Forcing|Anthropogenic|Albedo Change','Effective Radiative Forcing|Natural']):
    
    # returns the sum of specified driving rfs (by default those not included in FaIR):
    _out = RCMIP_forc.loc[('World',scenario,drivers)].droplevel(level=[0,1]).iloc[:,4:]
    _out = pd.DataFrame(_out.sum(axis=0,skipna=False).values,index=_out.columns.astype(int),columns=['forcing'])
    return _out

def RCMIP_to_FaIR_input_concs(scenario):
    _out = RCMIP_concs.loc[('World',scenario,FaIR_to_RCMIP_map.RCMIP_concs_key)].droplevel(level=[0,1]).iloc[:,4:]
    _out_index = _out.index.map({v: k for k, v in FaIR_to_RCMIP_map.RCMIP_concs_key.to_dict().items()})
    _out.set_index(_out_index,inplace=True)
    _out = _out.T.mul(FaIR_to_RCMIP_map.RCMIP_concs_scaling)
    _out.index = _out.index.astype(int)
    return _out.apply(pd.to_numeric)
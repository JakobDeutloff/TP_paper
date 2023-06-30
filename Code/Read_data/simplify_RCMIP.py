"""
Write out SSP scenario emmissions and forcings to save memory when reading the files
"""

from Code.FAIR2.fair.scripts.data_retrieval import *
import os

os.chdir(r'C:\Users\jakob\OneDrive\Dokumente\Dokumente Uni\Master\4.Semester\Research Project\ModelFramework')


# %%
# Function from Leach 2021 notebooks
def get_SSP_emms(ssp):
    emms = RCMIP_to_FaIR_input_emms(ssp).interpolate()
    rebase_species = ['so2', 'nox', 'co', 'nmvoc', 'bc', 'nh3', 'oc', 'nox_avi', 'methyl_bromide', 'methyl_chloride',
                      'chcl3', 'ch2cl2']
    emms.loc[:, rebase_species] -= emms.loc[1750, rebase_species]
    return emms


# Theses are only the SSP scenarios, there are also RCP and idealised runs available
choose_ssps = ['ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp370-lowNTCF-aerchemmip', 'ssp370-lowNTCF-gidden', 'ssp434',
               'ssp460', 'ssp534-over', 'ssp585']

# emms
SSP_emms_all = pd.concat([get_SSP_emms(x) for x in choose_ssps], axis=1, keys=choose_ssps)
SSP_emms_all.columns.names = ['Scenario', 'Variable']
SSP_emms_all.index.names = ['Year']
## forcing
SSP_forc_all = pd.concat([get_RCMIP_forc(x) for x in choose_ssps], axis=1, keys=choose_ssps)
SSP_forc_all.columns.names = ['Scenario', 'Variable']
SSP_forc_all.index.names = ['Year']

# %% Save to csv files
SSP_emms_all.to_csv('Data/RCMIP_simplified/SSP_emms.csv')
SSP_forc_all.to_csv('Data/RCMIP_simplified/SSP_forc.csv')

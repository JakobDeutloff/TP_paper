"""
Script to explore errors in CTEM calibration - Ts for PFTP and AMAZ (Fig. S3) and Cum. C Emissions for PFAT (Fig. S4)
in 2100
"""
import numpy as np
from Code.Read_data.read_simplified_RCMIP import *
import matplotlib.pyplot as plt
from Code.Calibrate.Functions import *
from Code.Read_data.Constants_new import *
from Code.FAIR2.fair.fair_runnter_ctem import run_FaIR_intanal, get_gas_parameter_defaults, \
    get_thermal_parameter_defaults


# %% Define functions

def load_forc(ssp):
    """
    Load forcing and parameters to run fair
    :param ssp: ssp
    :return: thermal_params_ens, gas_params_ens, emms, forc
    """
    ens_size = 10 * 10
    N_sets = ens_size
    # members
    members = np.arange(N_sets)
    members = [str(x) for x in members]

    # get emissions and forcing
    emms = SSP_emms[[ssp]]
    emms.columns = emms.columns.remove_unused_levels().droplevel()
    emms = pd.concat([emms] * N_sets, axis=1, keys=members)

    forc = SSP_forc[[ssp]]
    forc.columns = forc.columns.remove_unused_levels().droplevel()
    forc = pd.concat([forc] * N_sets, axis=1, keys=members)

    # Get gas and thermal params
    gas_params = get_gas_parameter_defaults()
    idx = gas_params.index
    columns = pd.MultiIndex.from_product([members, list(gas_params.columns.levels[1])])
    gas_params_ens = pd.DataFrame(index=idx, columns=columns, data=np.tile(gas_params.values, N_sets))

    thermal_params = get_thermal_parameter_defaults()
    idx = thermal_params.index
    columns = pd.MultiIndex.from_product([members, list(thermal_params.columns.levels[1])])
    thermal_params_ens = pd.DataFrame(index=idx, columns=columns, data=np.tile(thermal_params.values, N_sets))

    return thermal_params_ens, gas_params_ens, emms, forc


def plot_errors_K(ax, errors, P_3d, F_3d, maxi):
    """
    Plot errors in K for PFAT in year 2100
    :param ax: ax
    :param errors: 2D array of calculated errors in %
    :param P_3d: 2D array of thresholds
    :param F_3d: 2D array of F100 feedbacks
    :return: pcolormesh for colorbar
    """
    pm = ax.pcolormesh(P_3d, F_3d, errors, vmin=-maxi, vmax=0, cmap='hot')

    xticks = np.zeros(len(P_3d[0, :10]) - 1)
    for i in range(len(P_3d[0, :10]) - 1):
        xticks[i] = (P_3d[0, :10][i] + P_3d[0, :10][i + 1]) / 2

    yticks = np.zeros(len(F_3d[:10, 0]) - 1)
    for i in range(len(F_3d[:10, 0]) - 1):
        yticks[i] = (F_3d[:10, 0][i] + F_3d[:10, 0][i + 1]) / 2

    ax.set_xticks(xticks, minor=True)
    ax.set_yticks(yticks, minor=True)

    ax.xaxis.grid(True, which='minor', color='grey')
    ax.yaxis.grid(True, which='minor', color='grey')
    ax.set_yticks(F_3d[:10, 0], minor=False)
    ax.set_xticks(P_3d[0, :10], minor=False)
    ax.set_xticklabels(P_3d[0, :10], rotation=45)

    ax.set_xlabel('Threshold [°C]')
    ax.set_ylabel(r'$\mathrm{F_{100} ~ [GtC ~ yr^{-1}]}$')

    return pm


def calculate_error_K(result):
    """
    Calculate the relative error in cum. carbon emissions for PFAT in 2100 and 2300
    :param result: result_PFAT from Fair
    :return: 2D relative errors and 1D true values for 2100 and 2300
    """

    # Initialize arrays of right shape
    error100 = np.zeros((10, 10))
    error300 = np.zeros((10, 10))

    # find correct values
    true_val_100 = result['T'].loc[2100].values * TE_params_PFAT['F100']
    # set to max impact if reached
    true_val_100[true_val_100 > TE_params_PFAT['K'][:, 2]] = \
        TE_params_PFAT['K'][:, 2][true_val_100 > TE_params_PFAT['K'][:, 2]]
    # set to nan if not tipped
    tip_100 = result['tipping_elements'].xs(('PFAT', 'bool_tip'), level=(1, 2), axis=1).loc[2100].values
    true_val_100[tip_100 < 1] = np.nan

    # same for 2300
    true_val_300 = result['T'].loc[2300].values * TE_params_PFAT['F300']
    tip_300 = result['tipping_elements'].xs(('PFAT', 'bool_tip'), level=(1, 2), axis=1).loc[2300].values
    true_val_300[true_val_300 > TE_params_PFAT['K'][:, 2]] = \
        TE_params_PFAT['K'][:, 2][true_val_300 > TE_params_PFAT['K'][:, 2]]
    true_val_300[tip_300 < 1] = np.nan

    # extract calculated values
    calc_val_100 = result['tipping_elements'].xs(('PFAT', 'C_stock'), level=(1, 2), axis=1).loc[2100].values
    calc_val_300 = result['tipping_elements'].xs(('PFAT', 'C_stock'), level=(1, 2), axis=1).loc[2300].values

    error100.ravel()[:] = ((calc_val_100 - true_val_100) / true_val_100) * 100
    error300.ravel()[:] = ((calc_val_300 - true_val_300) / true_val_300) * 100

    return error100, error300, true_val_100, true_val_300


#  calculate errors of TS for PFTP and AMAZ
def calculate_errors(result, elem, Ts):
    """
    Calculates errors in Ts for PFATP and AMAZ
    """
    times = np.zeros((10, 10))

    for i in range(100):
        trg_year = find_trg_year(result['tipping_elements'].xs((str(i), elem, 'bool_tip'), level=(0, 1, 2), axis=1))
        thr_year = find_threshold(result['tipping_elements'].xs((str(i), elem, 'C_stock'), level=(0, 1, 2), axis=1),
                                  0.995)
        if thr_year is None or trg_year is None:
            times.ravel()[i] = np.nan
        else:
            times.ravel()[i] = thr_year - trg_year

    errors = ((times - Ts) / Ts) * 100  # error in %

    return errors


def plot_errors(ax, errors, P_3d, K_3d, maxi):
    """
    Plots errors in Ts for PFTP and AMAZ
    """

    pm = ax.pcolormesh(P_3d, K_3d, errors, vmin=-maxi, vmax=maxi, cmap='RdBu')

    xticks = np.zeros(len(P_3d[0, :10]) - 1)
    for i in range(len(P_3d[0, :10]) - 1):
        xticks[i] = (P_3d[0, :10][i] + P_3d[0, :10][i + 1]) / 2

    yticks = np.zeros(len(K_3d[:10, 0]) - 1)
    for i in range(len(K_3d[:10, 0]) - 1):
        yticks[i] = (K_3d[:10, 0][i] + K_3d[:10, 0][i + 1]) / 2

    ax.set_xticks(xticks, minor=True)
    ax.set_yticks(yticks, minor=True)

    ax.xaxis.grid(True, which='minor', color='grey')
    ax.yaxis.grid(True, which='minor', color='grey')
    ax.set_yticks(K_3d[:10, 0], minor=False)
    ax.set_xticks(P_3d[0, :10], minor=False)
    ax.set_xticklabels(P_3d[0, :10], rotation=45)

    return pm


# Run whole analysis for PFAT
def test_K_PFAT(ssp, TE_params_PFAT):
    """
    Whole routine using above functions to generate plot of rlative error in cum. emissions in 2100
    """
    thermal_params_ens, gas_params_ens, emms, forc = load_forc(ssp)

    result_PFAT = run_FaIR_intanal(emissions_in=emms,
                                   forcing_in=forc,
                                   thermal_parameters=thermal_params_ens,
                                   gas_parameters=gas_params_ens,
                                   show_run_info=False,
                                   TE_params=TE_params_PFAT)

    # Calculate errors
    error_100, error_300, true_val_100, true_val_300 = calculate_error_K(result_PFAT)

    # Plot 2100 error
    fig, ax = plt.subplots(1, 1, figsize=(5, 4), sharey='row')
    maxi = np.nanmax(np.abs(error_100))
    max_round = np.round(maxi, -1)
    pm_100 = plot_errors_K(ax, error_100, P_PFAT_ext, F_PFAT_100_exp, max_round)
    # Colorbar
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.12, 0.05, 0.75])
    cbar = fig.colorbar(pm_100, cax=cbar_ax, orientation='vertical', pad=0.01, label='Cumulative Emissions Error [%]')
    cbar.set_ticks(np.arange(-max_round, 0+10, 10))
    plt.savefig('Plots/Calibration/K_errors_100' + ssp + '.png', bbox_inches='tight', dpi=300)
    plt.show()

    # Plot 2300 error
    fig, ax = plt.subplots(1, 1, figsize=(5, 4), sharey='row')
    maxi = np.nanmax(np.abs(error_300))
    max_round = np.round(maxi, -1)
    pm_300 = plot_errors_K(ax, error_300, P_PFAT_ext, F_PFAT_300_exp, max_round)
    # Colorbar
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.12, 0.05, 0.75])
    cbar = fig.colorbar(pm_300, cax=cbar_ax, orientation='vertical', pad=0.01, label='Cumulative Emissions Error [%]')
    cbar.set_ticks(np.arange(-max_round, 0+10, 10))
    plt.savefig('Plots/Calibration/K_errors_300' + ssp + '.png', bbox_inches='tight', dpi=300)
    plt.show()


def test_TS_PFTP_AMAZ(ssp, TE_params_PFTP_min, TE_params_PFTP_mean, TE_params_PFTP_max, TE_params_AMAZ_min,
                      TE_params_AMAZ_mean, TE_params_AMAZ_max):
    """
    Routine to plot relative error in Ts for PFTP and AMAZ
    """

    thermal_params_ens, gas_params_ens, emms, forc = load_forc(ssp)
    # PFTP
    result_PFTP_min = run_FaIR_intanal(emissions_in=emms,
                                       forcing_in=forc,
                                       thermal_parameters=thermal_params_ens,
                                       gas_parameters=gas_params_ens,
                                       show_run_info=False,
                                       TE_params=TE_params_PFTP_min)

    result_PFTP_mean = run_FaIR_intanal(emissions_in=emms,
                                        forcing_in=forc,
                                        thermal_parameters=thermal_params_ens,
                                        gas_parameters=gas_params_ens,
                                        show_run_info=False,
                                        TE_params=TE_params_PFTP_mean)

    result_PFTP_max = run_FaIR_intanal(emissions_in=emms,
                                       forcing_in=forc,
                                       thermal_parameters=thermal_params_ens,
                                       gas_parameters=gas_params_ens,
                                       show_run_info=False,
                                       TE_params=TE_params_PFTP_max)

    # AMAZ
    result_AMAZ_min = run_FaIR_intanal(emissions_in=emms,
                                       forcing_in=forc,
                                       thermal_parameters=thermal_params_ens,
                                       gas_parameters=gas_params_ens,
                                       show_run_info=False,
                                       TE_params=TE_params_AMAZ_min)

    result_AMAZ_mean = run_FaIR_intanal(emissions_in=emms,
                                        forcing_in=forc,
                                        thermal_parameters=thermal_params_ens,
                                        gas_parameters=gas_params_ens,
                                        show_run_info=False,
                                        TE_params=TE_params_AMAZ_mean)

    result_AMAZ_max = run_FaIR_intanal(emissions_in=emms,
                                       forcing_in=forc,
                                       thermal_parameters=thermal_params_ens,
                                       gas_parameters=gas_params_ens,
                                       show_run_info=False,
                                       TE_params=TE_params_AMAZ_max)

    # calculate errors
    errors_PFTP_min = calculate_errors(result_PFTP_min, 'PFTP', Ts_PFTP[0])
    errors_PFTP_mean = calculate_errors(result_PFTP_mean, 'PFTP', Ts_PFTP[1])
    errors_PFTP_max = calculate_errors(result_PFTP_max, 'PFTP', Ts_PFTP[2])

    errors_AMAZ_min = calculate_errors(result_AMAZ_min, 'AMAZ', Ts_AMAZ[0])
    errors_AMAZ_mean = calculate_errors(result_AMAZ_mean, 'AMAZ', Ts_AMAZ[1])
    errors_AMAZ_max = calculate_errors(result_AMAZ_max, 'AMAZ', Ts_AMAZ[2])

    # Plot
    fig, axes = plt.subplots(2, 3, figsize=(10, 7), sharey='row')

    errors_PFTP = [errors_PFTP_min, errors_PFTP_mean, errors_PFTP_max]
    errors_AMAZ = [errors_AMAZ_min, errors_AMAZ_mean, errors_AMAZ_max]

    maxi = np.nanmax(np.abs([errors_AMAZ, errors_PFTP]))
    max_round = np.round(maxi, -1)
    for i in range(3):
        # PFTP
        pm_PFTP = plot_errors(axes[0, i], errors_PFTP[i], P_PFTP_ext, K_PFTP_ext, max_round)
        # AMAZ
        pm_AMAZ = plot_errors(axes[1, i], errors_AMAZ[i], P_AMAZ_ext, K_AMAZ_ext, max_round)

    # Colorbar
    fig.subplots_adjust(bottom=0.23)
    cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.02])
    cbar = fig.colorbar(pm_AMAZ, cax=cbar_ax, orientation='horizontal', pad=0.01, label='Timescale Error [%]',
                        extend='max')
    cbar.set_ticks(np.arange(-max_round, max_round+10, 10))

    # Set labels
    titles = [r'Min $H$', r'Mean $H$', r'Max $H$']
    for i in range(3):
        axes[0, i].set_title(titles[i])
        axes[1, i].set_xlabel('Threshold [°C]')
    axes[0, 0].set_ylabel('PFTP Impact [GtC]')
    axes[1, 0].set_ylabel('AMAZ Impact [GtC]')

    plt.savefig('Plots/Calibration/TS_errors_' + ssp + '.png', bbox_inches='tight', dpi=300)

    plt.show()


# %% PFTP setup - Check error in Ts

ens_size = 100

# tipping threshold in °C
P_PFTP = np.linspace(P_min['PFTP'], P_max['PFTP'], 10).round(2)
P_AMAZ = np.array(ens_size * [100])  # disable AMAZ
P_PFAT = np.array(ens_size * [100])  # disable PFAT

# timescale
Ts_PFTP = [Ts_min['PFTP'], Ts_mean['PFTP'], Ts_max['PFTP']]

# Rate
R_PFTP = np.interp(Ts_PFTP, r_PFTP.timescale, r_PFTP.rate).squeeze()
R_AMAZ = np.array(ens_size * [0.1615])
R_PFAT = np.array(ens_size * [0.1615])

# impact in GtC
K_PFTP = np.linspace(K_min['PFTP'], K_max['PFTP'], 10).round()
K_AMAZ = np.array(ens_size * [0])
K_PFAT = np.array(ens_size * [0])

# Create meshgrid to cover all possible variations
P_PFTP_ext, K_PFTP_ext = np.meshgrid(P_PFTP, K_PFTP)
P_PFTP = P_PFTP_ext.reshape(100)
K_PFTP = K_PFTP_ext.reshape(100)

# Initial stock
y0_PFTP = 0.005 * K_PFTP
y0_AMAZ = np.array(ens_size * [0])
y0_PFAT = np.array(ens_size * [0])

# Feedback of PFAT - not used here because we calibrate it with mean impact
F_PFAT_100 = np.array(ens_size * [0])
F_PFAT_300 = np.array(ens_size * [0])
cal_flag = False  # only used for calibration of PFAT

# arange parameter sets for all three timescales
TE_params_PFTP_min = {'P': np.array([P_PFTP, P_AMAZ, P_PFAT]).T, 'R': np.array([100 * [R_PFTP[0]], R_AMAZ, R_PFAT]).T,
                      'K': np.array([K_PFTP, K_AMAZ, K_PFAT]).T,
                      'CH4': np.array([CH4_frac['PFTP'], 0, CH4_frac['PFAT']]),
                      'y0': np.array([y0_PFTP, y0_AMAZ, y0_PFAT]), 'F100': F_PFAT_100, 'F300': F_PFAT_300,
                      'cal_flag': cal_flag}

TE_params_PFTP_mean = {'P': np.array([P_PFTP, P_AMAZ, P_PFAT]).T, 'R': np.array([100 * [R_PFTP[1]], R_AMAZ, R_PFAT]).T,
                       'K': np.array([K_PFTP, K_AMAZ, K_PFAT]).T,
                       'CH4': np.array([CH4_frac['PFTP'], 0, CH4_frac['PFAT']]),
                       'y0': np.array([y0_PFTP, y0_AMAZ, y0_PFAT]), 'F100': F_PFAT_100, 'F300': F_PFAT_300,
                       'cal_flag': cal_flag}

TE_params_PFTP_max = {'P': np.array([P_PFTP, P_AMAZ, P_PFAT]).T, 'R': np.array([100 * [R_PFTP[2]], R_AMAZ, R_PFAT]).T,
                      'K': np.array([K_PFTP, K_AMAZ, K_PFAT]).T,
                      'CH4': np.array([CH4_frac['PFTP'], 0, CH4_frac['PFAT']]),
                      'y0': np.array([y0_PFTP, y0_AMAZ, y0_PFAT]), 'F100': F_PFAT_100, 'F300': F_PFAT_300,
                      'cal_flag': cal_flag}

# %% AMAZ setup - Check error in Ts

# tipping threshold in °C
P_AMAZ = np.linspace(P_min['AMAZ'], P_max['AMAZ'], 10).round(2)
P_PFTP = np.array(ens_size * [100])  # disable PFTP
P_PFAT = np.array(ens_size * [100])  # disable PFAT

# timescale
Ts_AMAZ = [Ts_min['AMAZ'], Ts_mean['AMAZ'], Ts_max['AMAZ']]

# Rate
R_AMAZ = np.interp(Ts_AMAZ, r_AMAZ.timescale, r_AMAZ.rate).squeeze()
R_PFTP = np.array(ens_size * [0.1615])
R_PFAT = np.array(ens_size * [0.1615])

# impact in GtC
K_AMAZ = np.linspace(K_min['AMAZ'], K_max['AMAZ'], 10).round()
K_PFTP = np.array(ens_size * [0])
K_PFAT = np.array(ens_size * [0])

# Create meshgrid to cover all possible variations
P_AMAZ_ext, K_AMAZ_ext = np.meshgrid(P_AMAZ, K_AMAZ)
P_AMAZ = P_AMAZ_ext.reshape(100)
K_AMAZ = K_AMAZ_ext.reshape(100)

# Initial stock
y0_AMAZ = 0.005 * K_AMAZ
y0_PFTP = np.array(ens_size * [0])
y0_PFAT = np.array(ens_size * [0])

# Feedback of PFAT - not used here because we calibrate it with mean impact
F_PFAT_100 = np.array(ens_size * [0])
F_PFAT_300 = np.array(ens_size * [0])
cal_flag = False  # only used for calibration of PFAT

# arange parameter sets for all three timescales
TE_params_AMAZ_min = {'P': np.array([P_PFTP, P_AMAZ, P_PFAT]).T, 'R': np.array([R_PFTP, 100 * [R_AMAZ[0]], R_PFAT]).T,
                      'K': np.array([K_PFTP, K_AMAZ, K_PFAT]).T,
                      'CH4': np.array([CH4_frac['PFTP'], 0, CH4_frac['PFAT']]),
                      'y0': np.array([y0_PFTP, y0_AMAZ, y0_PFAT]), 'F100': F_PFAT_100, 'F300': F_PFAT_300,
                      'cal_flag': cal_flag}

TE_params_AMAZ_mean = {'P': np.array([P_PFTP, P_AMAZ, P_PFAT]).T, 'R': np.array([R_PFTP, 100 * [R_AMAZ[1]], R_PFAT]).T,
                       'K': np.array([K_PFTP, K_AMAZ, K_PFAT]).T,
                       'CH4': np.array([CH4_frac['PFTP'], 0, CH4_frac['PFAT']]),
                       'y0': np.array([y0_PFTP, y0_AMAZ, y0_PFAT]), 'F100': F_PFAT_100, 'F300': F_PFAT_300,
                       'cal_flag': cal_flag}

TE_params_AMAZ_max = {'P': np.array([P_PFTP, P_AMAZ, P_PFAT]).T, 'R': np.array([R_PFTP, 100 * [R_AMAZ[2]], R_PFAT]).T,
                      'K': np.array([K_PFTP, K_AMAZ, K_PFAT]).T,
                      'CH4': np.array([CH4_frac['PFTP'], 0, CH4_frac['PFAT']]),
                      'y0': np.array([y0_PFTP, y0_AMAZ, y0_PFAT]), 'F100': F_PFAT_100, 'F300': F_PFAT_300,
                      'cal_flag': cal_flag}

# %% PFAT setup - check error in K

# tipping threshold in °C
P_PFTP = np.array(ens_size * [100])  # disable PFTP
P_AMAZ = np.array(ens_size * [100])  # disable AMAZ
P_PFAT = np.linspace(P_min['PFAT'], P_max['PFAT'], 10).round(2)

# timescale
Ts_PFAT = Ts_mean['PFAT']

# Rate
R_PFTP = np.array(ens_size * [0.1615])
R_AMAZ = np.array(ens_size * [0.1615])
R_PFAT = np.array(ens_size * [r_PFAT.rate[0]])

# impact in GtC
K_PFTP = np.array(ens_size * [0])
K_AMAZ = np.array(ens_size * [0])
K_PFAT = np.linspace(K_min['PFAT'], K_max['PFAT'], 10).round()

# Feedback of PFAT - corresponds to K
F_PFAT_100 = np.linspace(F_min['PFAT100'], F_max['PFAT100'], 10).round(1)
F_PFAT_300 = np.linspace(F_min['PFAT300'], F_max['PFAT300'], 10).round(1)
cal_flag = False

# Create meshgrid to cover all possible variations
P_PFAT_ext, K_PFAT_ext = np.meshgrid(P_PFAT, K_PFAT)
a, F_PFAT_100_exp = np.meshgrid(P_PFAT, F_PFAT_100)
a, F_PFAT_300_exp = np.meshgrid(P_PFAT, F_PFAT_300)
P_PFAT = P_PFAT_ext.reshape(100)
K_PFAT = K_PFAT_ext.reshape(100)
F_PFAT_100 = F_PFAT_100_exp.reshape(100)
F_PFAT_300 = F_PFAT_300_exp.reshape(100)

# Initial stock
y0_PFTP = np.array(ens_size * [0])
y0_AMAZ = np.array(ens_size * [0])
y0_PFAT = 0.005 * K_PFAT

# arange parameter sets
TE_params_PFAT = {'P': np.array([P_PFTP, P_AMAZ, P_PFAT]).T, 'R': np.array([R_PFTP, R_AMAZ, R_PFAT]).T,
                  'K': np.array([K_PFTP, K_AMAZ, K_PFAT]).T,
                  'CH4': np.array([CH4_frac['PFTP'], 0, CH4_frac['PFAT']]),
                  'y0': np.array([y0_PFTP, y0_AMAZ, y0_PFAT]), 'F100': F_PFAT_100, 'F300': F_PFAT_300,
                  'cal_flag': cal_flag}

# %% run analysis for selected ssps

ssps = ['ssp585']

for ssp in ssps:
    test_K_PFAT(ssp, TE_params_PFAT)
    #test_TS_PFTP_AMAZ(ssp, TE_params_PFTP_min, TE_params_PFTP_mean, TE_params_PFTP_max, TE_params_AMAZ_min,
    #                  TE_params_AMAZ_mean, TE_params_AMAZ_max)

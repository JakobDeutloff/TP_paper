"""
Old plotscripts used in the climate study project
"""
import matplotlib.pyplot as plt
import matplotlib
font = {'size': 10}
matplotlib.rc('font', **font)

def plot_fair_int_std(run_int, run_std, out_TE):

    """
    Plot to compare differences caused by constant alpha
    """

    col_int = 'red'
    col_std = 'blue'

    fig, axis = plt.subplots(2, 3, sharex=True, figsize=(14, 7))

    # Plot run with TE
    run_int['T'].plot(ax=axis[0, 0], legend=False, color=col_int)
    axis[0, 0].set_ylabel('T [°C]')
    run_int['C'].xs('carbon_dioxide', axis=1, level=-1).plot(ax=axis[0, 1], legend=False, color=col_int)
    axis[0, 1].set_ylabel('CO2_conc [ppm]')
    run_int['C'].xs('methane', axis=1, level=-1).plot(ax=axis[0, 2], legend=False, color=col_int)
    axis[0, 2].set_ylabel('CH4_conc [ppt]')
    run_int['alpha'].xs('carbon_dioxide', axis=1, level=-1).plot(ax=axis[1, 1], legend=False, color=col_int)
    axis[1, 0].set_ylabel('CO2_emm []')
    run_int['alpha'].xs('methane', axis=1, level=-1).plot(ax=axis[1, 2], legend=False, color=col_int)
    axis[1, 1].set_ylabel('alpha_CO2')
    run_int['Emissions'].xs('carbon_dioxide', axis=1, level=-1).plot(ax=axis[1, 0], legend=False, color=col_std)
    axis[1, 0].plot(run_int['T'].index, out_TE['emmAMAZ'], color='g')
    axis[1, 2].set_ylabel('alpha_CH4')

    # Plot run without TE
    run_std['T'].plot(ax=axis[0, 0], legend=False, color=col_std, linestyle='--')
    axis[0, 0].set_ylabel('T [°C]')
    run_std['C'].xs('carbon_dioxide', axis=1, level=-1).plot(ax=axis[0, 1], legend=False, color=col_std, linestyle='--')
    axis[0, 1].set_ylabel('CO2_conc [ppm]')
    run_std['C'].xs('methane', axis=1, level=-1).plot(ax=axis[0, 2], legend=False, color=col_std, linestyle='--')
    axis[0, 2].set_ylabel('CH4_conc [ppt]')
    run_std['alpha'].xs('carbon_dioxide', axis=1, level=-1).plot(ax=axis[1, 1], legend=False, color=col_std,
                                                                 linestyle='--')
    axis[1, 0].set_ylabel('CO2_emm []')
    run_std['alpha'].xs('methane', axis=1, level=-1).plot(ax=axis[1, 2], legend=False, color=col_std, linestyle='--')
    axis[1, 1].set_ylabel('alpha_CO2')
    run_std['Emissions'].xs('carbon_dioxide', axis=1, level=-1).plot(ax=axis[1, 0], legend=False, color=col_std,
                                                                     linestyle='--')
    axis[1, 2].set_ylabel('alpha_CH4')

    for ax in axis.flat:
        ax.grid()

    fig.suptitle('Scenario: ' + run_int['C'].columns[0][0] + ', Parameters: ' + run_int['C'].columns[0][1])

    plt.tight_layout()
    return fig


def plot_output_all_TE(run_int, run_std, out_TE, title):

    """
    Plot differences between standart and interactive run in therms of T, emmC, C and stocks ot TEs
    """

    col_int = 'red'
    col_std = 'blue'

    fig, axes = plt.subplots(2, 2, sharex='col', figsize=(10, 7))

    # Temp
    axes[0, 0].plot(run_int['T'], color=col_int, label='Interactive')
    axes[0, 0].plot(run_std['T'], color=col_std, label='Standard')
    axes[0, 0].set_ylabel('T [°C]')

    # CO2 conc
    axes[0, 1].plot(run_int['C'].xs('carbon_dioxide', axis=1, level=-1), color=col_int, label='Interactive')
    axes[0, 1].plot(run_std['C'].xs('carbon_dioxide', axis=1, level=-1), color=col_std, label='Standard')
    axes[0, 1].set_ylabel('CO2_conc [ppm]')

    # Emissions CO2
    axes[1, 0].plot(run_int['Emissions'].xs('carbon_dioxide', axis=1, level=-1), color='k', label='SSP')
    axes[1, 0].plot(out_TE['emmPFTP_CO2'], color='darkgreen', linestyle='--', label='PFTP')
    axes[1, 0].plot(out_TE['emmPFAT_CO2'], color='darkgreen', linestyle=':', label='PFAT')
    axes[1, 0].plot(out_TE['emmPFGT_CO2'], color='darkgreen', linestyle='-.', label='PFGT')
    axes[1, 0].set_ylabel('Carbon Emissions [GtC/yr]')

    # Emissions CH4
    axes[1, 1].plot(run_int['Emissions'].xs('methane', axis=1, level=-1), color='k', label='SSP')
    axes[1, 1].plot(out_TE['emmPFTP_CH4'], color='darkgreen', linestyle='--', label='PFTP')
    axes[1, 1].plot(out_TE['emmPFAT_CH4'], color='darkgreen', linestyle=':', label='PFAT')
    axes[1, 1].plot(out_TE['emmPFGT_CH4'], color='darkgreen', linestyle='-.', label='PFGT')
    axes[1, 1].set_ylabel('Carbon Emissions [MtCH4/yr]')


    for ax in axes.flat:
        ax.legend()
        ax.grid()

    fig.suptitle(title)
    plt.tight_layout()

    return fig


def plot_output_all_TE_2(run_int, run_std, out_TE, title):

    """
    Plot differences between standart and interactive run in therms of T, emmCO2, emmCH4 and stocks ot TEs
    """

    col_int = 'blue'
    col_std = 'black'

    fig, axes = plt.subplots(2, 2, sharex='col', figsize=(8, 6))

    # Temp
    axes[0, 0].plot(run_int['T'], color=col_int, label='Permafrost')
    axes[0, 0].plot(run_std['T'], color=col_std, linestyle='--', label='Standard')
    axes[0, 0].set_ylabel('GMST [°C]')

    # Stocks
    axes[0, 1].plot(out_TE['sPFTP'], color='red', linestyle='--', label='PFTP')
    axes[0, 1].plot(out_TE['sPFAT'], color='orange', linestyle=':', label='PFAT')
    axes[0, 1].plot(out_TE['sPFGT'], color='green', linestyle='-.', label='PFGT')
    axes[0, 1].set_ylabel('Cumulative C Emissions [GtC]')

    # Emissions CO2
    axes[1, 0].plot(run_int['Emissions'].xs('carbon_dioxide', axis=1, level=-1), color='k', label='SSP')
    axes[1, 0].plot(out_TE['emmPFTP_CO2'], color='red', linestyle='--', label='PFTP')
    axes[1, 0].plot(out_TE['emmPFAT_CO2'], color='orange', linestyle=':', label='PFAT')
    axes[1, 0].plot(out_TE['emmPFGT_CO2'], color='green', linestyle='-.', label='PFGT')
    axes[1, 0].set_ylabel(r'$\mathrm{CO_2}$ Emissions [GtC/yr]')

    # Emissions CH4
    axes[1, 1].plot(run_int['Emissions'].xs('methane', axis=1, level=-1), color='k', label='SSP')
    axes[1, 1].plot(out_TE['emmPFTP_CH4'], color='red', linestyle='--', label='PFTP')
    axes[1, 1].plot(out_TE['emmPFAT_CH4'], color='orange', linestyle=':', label='PFAT')
    axes[1, 1].plot(out_TE['emmPFGT_CH4'], color='green', linestyle='-.', label='PFGT')
    axes[1, 1].set_ylabel(r'$\mathrm{CH_4}$ Emissions [MtCH4/yr]')


    for ax in axes.flat:
        ax.grid()

    # handles, labels = axes[0, 1].get_legend_handles_labels()
    handles, labels = [(a + b) for a, b in zip(axes[0,0].get_legend_handles_labels(), axes[1, 0].get_legend_handles_labels())]
    fig.legend(handles, labels, loc='lower center', borderaxespad=0.1, ncol=3)
    fig.suptitle(title)
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.13)


    return fig
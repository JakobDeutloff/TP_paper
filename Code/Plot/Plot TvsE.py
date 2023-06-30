from Code.Read_data.Read_SSP_output import read_rand_mems, read_SSP_outout
import matplotlib.pyplot as plt


# %% def plot function

def plot_TvsE(ssp):

    end_year = 2500
    # Calculate cum Emissions
    weight_C = 12.011  # weight of carbon in u
    weight_CH4 = 16  # weight of CH4 in u
    fac_CH4 = (weight_C / weight_CH4) * 1e3  # used in fair for conversion from GtC to MtCH4
    em_const = res_const['Emm'][ssp].loc[:end_year, 'carbon_dioxide'] + res_const['Emm'][ssp].loc[:end_year, 'methane'] * (
            1 / fac_CH4)
    cum_em_const = em_const.cumsum()

    fig, ax = plt.subplots(figsize=(5, 5))

    # Hist
    ax.plot(cum_em_const.loc[:2016], res_const['T'][ssp].loc[:2016, '0.5'], color='grey')
    ax.plot(cum_em_const.loc[:2016], res_const['T'][ssp].loc[:2016, '0.95'], color='grey', linewidth=0.5)
    ax.plot(cum_em_const.loc[:2016], res_const['T'][ssp].loc[:2016, '0.05'], color='grey', linewidth=0.5)
    ax.fill_between(cum_em_const.loc[:2016], res_const['T'][ssp].loc[:2016, '0.05'],
                    res_const['T'][ssp].loc[:2016, '0.95'],
                    color='grey', alpha=0.3)

    # Const
    ax.plot(cum_em_const.loc[2016:], res_const['T'][ssp].loc[2016:end_year, '0.5'], color='k')
    ax.plot(cum_em_const.loc[2016:], res_const['T'][ssp].loc[2016:end_year, '0.95'], color='k', linewidth=0.5)
    ax.plot(cum_em_const.loc[2016:], res_const['T'][ssp].loc[2016:end_year, '0.05'], color='k', linewidth=0.5)
    ax.fill_between(cum_em_const.loc[2016:], res_const['T'][ssp].loc[2016:end_year, '0.05'],
                    res_const['T'][ssp].loc[2016:end_year, '0.95'],
                    color='grey', alpha=0.3)

    # Tip
    ax.plot(cum_em_const.loc[2016:], res_tip['T'][ssp].loc[2016:end_year, '0.5'], color='blue')
    ax.plot(cum_em_const.loc[2016:], res_tip['T'][ssp].loc[2016:end_year, '0.95'], color='blue', linewidth=0.5)
    ax.plot(cum_em_const.loc[2016:], res_tip['T'][ssp].loc[2016:end_year, '0.05'], color='blue', linewidth=0.5)
    ax.fill_between(cum_em_const.loc[2016:], res_tip['T'][ssp].loc[2016:end_year, '0.05'],
                    res_tip['T'][ssp].loc[2016:end_year, '0.95'],
                    color='blue', alpha=0.3)

    # Rand Members
    ax.plot(cum_em_const.loc[2016:], rand_mems[ssp]['T'].loc[2016:end_year].iloc[:, [1, 2, 16, 17, 18, 20, 30, 31]],
            color='blue',
            linewidth=0.7, linestyle='--')

    ax.set_xlabel('Cumulative Emissions [GtC]')
    ax.set_ylabel('GMST [°C]')
    ax.set_title(ssp)
    plt.tight_layout()
    plt.savefig('Plots/Comparison/5000_tip_2_5000_const/TvsE/' + ssp + '_TvsE.png', dpi=500)
    plt.show()


if __name__ == '__main__':

    # Load results
    rand_mems = read_rand_mems('5000_tip_2')
    res_tip = read_SSP_outout('5000_tip_2')
    res_const = read_SSP_outout('5000_const')
    for ssp in list(res_tip['T'].columns.levels[0]):
        plot_TvsE(ssp)

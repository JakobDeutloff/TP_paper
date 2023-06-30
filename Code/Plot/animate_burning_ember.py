from Code.Read_data.Read_SSP_output import read_SSP_outout, read_probs
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import matplotlib as mpl
from matplotlib import animation

# %% Read data
ens_name_tip = '5000_tip_2'
ens_name_const = '5000_const_2'

SSP_out_tip = read_SSP_outout(ens_name_tip)
SSP_out_const = read_SSP_outout(ens_name_const)
Probs_tip = read_probs(ens_name_tip)
Probs_const = read_probs(ens_name_const)
# %% Define functions

def gradientbars(bars, cmap):
    ax = bars[0].axes
    lim = ax.get_xlim() + ax.get_ylim()
    ax.axis(lim)
    for bar in bars:
        bar.set_facecolor("none")
        x, y = bar.get_xy()
        w, h = bar.get_width(), bar.get_height()
        grad = np.atleast_2d(np.linspace(0, 1 * h / 1, 1000)).T
        # zorder of 2 to get gradients above the facecolor, but below the bar outlines
        ax.imshow(grad, extent=[x, x + w, y, y + h], origin='lower', aspect="auto", zorder=2,
                  norm=cm.colors.NoNorm(vmin=0, vmax=1), cmap=plt.get_cmap(cmap))

def plot_only_burning_ember(ssp, year, ax, columns):


    values_max_const = SSP_out_const['Tip_probs_total'].loc[year][ssp][columns].values
    values_max_tip = SSP_out_tip['Tip_probs_total'].loc[year][ssp][columns].values

    # set grid and label of ax
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.05, 0.33, 0.5, 0.66, 0.95, 1])
    ax.grid(which='major', axis='y', linestyle='--', color='gray', zorder=0)

    # Plot bars of coupled run
    ax.bar(columns, values_max_tip, color='mediumvioletred', edgecolor='grey', zorder=2)
    # plot bars of standard run
    my_bar = ax.bar(columns, values_max_const, color='grey', edgecolor='grey', zorder=3)
    # Make fancy BE optics
    gradientbars(my_bar, 'YlOrRd')
    # Title
    ax.set_title(ssp + ', Year: ' + str(year))

    # modify xticks
    xticks_loc = ax.get_xticks()
    ax.set_xticks(xticks_loc)
    ax.set_xticklabels(columns, rotation=0, fontdict={'fontsize': 7})
    ax.set_ylabel('Probability of Tipping')


def plot_T_ssp(ssp, ax, year):

    # set extent
    ax.set_xlim([2014, 2500])
    ax.set_ylim([0, T_max])

    # plot temperature uncoupled
    ax.plot(SSP_out_const['T'].loc[2014:year, (ssp, '0.5')], color='k', label=ssp)
    ax.fill_between(np.arange(2014, year+1), *SSP_out_const['T'].loc[2014:year, (ssp, ['0.05', '0.95'])].values.T,
                    color='k', alpha=0.3, linestyle='-')

    # plot temp and quantiles coupled
    ax.plot(SSP_out_tip['T'].loc[2014:year, (ssp, '0.5')], color='r', label=ssp)
    ax.fill_between(np.arange(2014, year+1), *SSP_out_tip['T'].loc[2014:year, (ssp, ['0.05', '0.95'])].values.T,
                    color='r', alpha=0.3, linestyle='-')


    # set labes
    ax.set_xlabel('Year')
    ax.set_ylabel('GMST anomaly [°C]')





# %% Run animation

# set path to ffmpeg programm
mpl.rcParams['animation.ffmpeg_path'] = \
    r'C:\Program Files\ffmpeg\bin\ffmpeg.exe'

years = np.arange(2015, 2501)
columns = Probs_const['Years'].columns.levels[1]


# define animation
def animate(year):
    # burning ember
    axes[0].clear()
    plot_only_burning_ember(ssp, year, axes[0], columns)
    # temperature
    axes[1].clear()
    plot_T_ssp(ssp, axes[1], year)


for ssp in ['ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585']:

    # create figure
    fig, axes = plt.subplots(2, 1, gridspec_kw={'height_ratios':(2, 1)}, figsize=(8, 7))
    # get max temp for ssp
    T_max = SSP_out_tip['T'][ssp]['0.95'].max()

    # create animation
    anim = animation.FuncAnimation(fig,
                                   animate,
                                   years,
                                   interval=100,
                                   blit=False, repeat=True)

    # save animation
    writer = animation.FFMpegWriter(fps=10)
    anim.save('Plots/SSP_ensembles/5000_tip_2/Animation/' + ssp + '.mp4', dpi=200, writer=writer)

"""
Create Burning Ember plots (Fig. 6, 7, 8, S11 - S17) and calculate how much earlier tipping happens on average due to carbon TEs
"""

# %%
from Code.Read_data.Read_SSP_output import read_SSP_outout, read_probs
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import pandas as pd

# %% Read data
ens_name_tip = "coupled_ensemble"
ens_name_const = "uncoupled_ensemble"

SSP_out_tip = read_SSP_outout(ens_name_tip)
Probs_tip = read_probs(ens_name_tip)
SSP_out_const = read_SSP_outout(ens_name_const)
Probs_const = read_probs(ens_name_const)


# %% Define Functions


def plot_burning_emeber(ssp, max_year=2500):
    columns = Probs_const["Years"].columns.levels[1]
    values_max_const = (
        SSP_out_const["Tip_probs_total"].loc[:max_year][ssp][columns].max().values
    )
    values_max_tip = (
        SSP_out_tip["Tip_probs_total"].loc[:max_year][ssp][columns].max().values
    )
    values_slow_tipping = Probs_const["slow_tip_probs"].loc[2450, ssp][columns].values

    fig, ax = plt.subplots(figsize=(8, 5))

    # set grid
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.1, 0.5, 0.9, 1])
    ax.grid(which="major", axis="y", linestyle="--", color="gray", zorder=0)

    # Plot bars of coupled run
    bars1 = ax.bar(
        columns,
        values_max_tip,
        color="#E50000",
        edgecolor="grey",
        zorder=2,
    )
    # plot instantaneous tipping probabilities (uncoupled)
    bars2 = ax.bar(
        columns,
        values_max_const,
        color="#FFA500",
        edgecolor="grey",
        zorder=3,
    )
    # plot slow tipping probabilities (uncoupled)
    if max_year == 2500:
        bars3 = ax.bar(
            columns,
            values_slow_tipping,
            color="#FFFF14",
            edgecolor="grey",
            zorder=4,
        )
    first, last = ax.get_xlim()
    plt.xlim(left=first + 0.7, right=last - 0.7)
    ax.set_ylabel("Probability of Triggering")
    # Create legend
    ax.legend(
        [bars3[0], bars2[0], bars1[0]],
        ["Delayed", "Instantaneous", "Instantaneous + Carbon TEs"],
        loc="lower center",
        bbox_to_anchor=(0.5, -0.17),
        ncol=3,
    )

    ax.set_title(ssp)
    plt.tight_layout()


# %% Plot Burning Embers diagramm

for ssp in Probs_const["Years"].columns.levels[0]:
    plot_burning_emeber(ssp, 2500)
    plt.tight_layout()
    plt.savefig("Plots/Burning_ember_2500/" + ssp + ".png", dpi=300)

# %% calculate total increase in probability of tippig for each SSP
ssps = list(Probs_const["Years"].columns.levels[0])

P_inst = pd.DataFrame(index=[0], columns=ssps)
P_inst_TEs = pd.DataFrame(index=[0], columns=ssps)
P_delayed = pd.DataFrame(index=[0], columns=ssps)
P_inst_200 = pd.DataFrame(index=[0], columns=ssps)

for ssp in ssps:
    P_inst.loc[:, (ssp)] = SSP_out_const["Tip_probs_total"][ssp].max().mean().round(2)
    P_inst_200.loc[:, (ssp)] = (
        SSP_out_const["Tip_probs_total"][ssp].loc[:2200].max().mean().round(2)
    )
    P_inst_TEs.loc[:, (ssp)] = SSP_out_tip["Tip_probs_total"][ssp].max().mean().round(2)
    P_delayed.loc[:, (ssp)] = (
        Probs_const["slow_tip_probs"].loc[2450, ssp].mean().round(2)
    )


# %%

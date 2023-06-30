"""
Quick script contining the numbers needed for CH4 and CO2 conversions and some calculations
"""

import numpy as np

GWP_CH4 = 27.9
m_CH4 = 16  # in u
m_CO2 = 44  # in u
m_C = 12  # in u


def CH4toCO2e(CH4):
    """
    takes CH4 in GtC to give CO2e in GtC
    :param CH4: CH4 emm in GtC
    :return: CO2e in GtCO2
    """
    # Get CO2e in GtC
    CO2e = CH4 * (m_CH4 / m_C) * GWP_CH4

    return CO2e


def CH4toC(CH4):
    """
    :param CH4: CH4 in TgCH4
    :return: C in GtC
    """
    return CH4 * 1e-3 * (m_C / m_CH4)


def CH4CeqtoC(CH4):
    """
    :param CH4: CH4 in GtCeq
    :return: CH4 in GtC
    """
    return CH4 * (m_CO2 / (GWP_CH4 * m_CH4))


# %% Values Turetsky
CH4_100_85 = 1.180  # in PgC - GtC
CH4_300_85 = 1.970
CO2_100_85 = 3.1  # in GtC
CO2_300_85 = 7.2

CH4_100_45 = 2.330
CH4_300_45 = 5.605
CO2_100_45 = 2.3
CO2_300_45 = 11.6

# %% calculate average CH4 percentage

CH4 = np.array([CH4_100_85, CH4_300_85, CH4_100_45, CH4_300_45])
CO2 = np.array([CO2_100_85, CO2_300_85, CO2_100_45, CO2_300_45])

perc = np.mean(CH4 / (CO2)) * 100
print(perc)

# %% calculate amplifcations

ampli = CH4toCO2e(2.3) + 97.7 / 100
print(ampli)

# %% convert Schneider von Deimling 2015 values
# CO2 valus in GtC for 2100, 2200, 2300 with uncertainty ranges (first value mean second and third lower,
# upper uncertainty range)
RCP26_CO2 = [36, 20, 58, 56, 35, 89, 64, 40, 98]
RCP45_CO2 = [54, 28, 92, 118, 75, 180, 155, 104, 216]
RCP60_CO2 = [60, 29, 101, 156, 103, 224, 193, 134, 270]
RCP85_CO2 = [87, 42, 141, 194, 136, 270, 228, 157, 313]

# CH4 values in TgCH4
RCP26_CH4 = [446, 218, 921, 818, 410, 1753, 1035, 539, 2236]
RCP45_CH4 = [1126, 538, 2356, 3117, 1657, 5969, 4705, 2592, 8449]
RCP60_CH4 = [1270, 663, 2440, 3104, 1818, 5372, 4615, 2664, 7778]
RCP85_CH4 = [1474, 836, 2614, 3592, 2141, 6093, 5877, 3644, 9989]

# C values from Gasser
RCP26_C_G = [27, 39, 47]
RCP45_C_G = [35, 64, 89]
RCP60_C_G = [42, 99, 145]
RCP85_C_G = [59, 150, 212]

# %% Calculate GtC for Schneider
RCP26_C = RCP26_CO2 + CH4toC(np.array(RCP26_CH4))
RCP45_C = RCP45_CO2 + CH4toC(np.array(RCP45_CH4))
RCP60_C = RCP60_CO2 + CH4toC(np.array(RCP60_CH4))
RCP85_C = RCP85_CO2 + CH4toC(np.array(RCP85_CH4))

# %% compare estimates
perc_26 = RCP26_C[[0, 3, 6]] / RCP26_C_G
perc_45 = RCP45_C[[0, 3, 6]] / RCP45_C_G
perc_60 = RCP60_C[[0, 3, 6]] / RCP60_C_G
perc_85 = RCP85_C[[0, 3, 6]] / RCP85_C_G
perc_total = np.mean([perc_26, perc_45, perc_60, perc_85])

# %% translate IPCC feedback to carbon

IPCC_CO2_min = 3.1  # in GtC
IPCC_CO2_mean = 18
IPCC_CO2_max = 41
IPCC_CH4_min = 0.7  # in GtCeq
IPCC_CH4_mean = 2.8
IPCC_CH4_max = 7.3

IPCC_CH4_min_C = CH4CeqtoC(IPCC_CH4_min)
IPCC_CH4_mean_C = CH4CeqtoC(IPCC_CH4_mean)
IPCC_CH4_max_C = CH4CeqtoC(IPCC_CH4_max)

C_min = IPCC_CH4_min_C + IPCC_CO2_min
C_mean = IPCC_CH4_mean_C + IPCC_CO2_mean
C_max = IPCC_CH4_max_C + IPCC_CO2_max

CH4_frac = (IPCC_CH4_mean_C / C_mean) * 100

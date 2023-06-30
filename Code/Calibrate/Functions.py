"""
Contains the functions used in the calibration scripts
"""


# Find trigger time:
def find_trg_year(trg):
    for i in range(len(trg)):
        if trg.iloc[i].values[0] > 0:
            return trg.index[i]


# find threshold time
def find_threshold(S, thr):
    for i in range(len(S)):
        if S.iloc[i].values[0] > thr * S.max().values[0]:
            return S.index[i]


# find year when T_goal is exceeded
def find_temp(T, T_goal):
    for i in range(len(T)):
        if T.iloc[i].values[0] > T_goal:
            year_tip = T.index[i]
            break
    return year_tip

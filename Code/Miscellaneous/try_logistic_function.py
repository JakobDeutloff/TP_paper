from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt
import os
os.chdir(r'/')


def model(t, y, R, Im):
    dy = R * y * (1 - y / Im)

    return dy


R = 10/100
Im = 100
y0 = [0.01]
out = integrate.solve_ivp(model, y0=y0, t_span=(0, 120), t_eval=np.arange(0, 120), args=[R, Im])

# plot output
fig, ax = plt.subplots()
ax.plot(out.t, out.y.squeeze())
ax.axvline(100, color='red', linestyle='--', label='Timescale')
ax.legend()
ax.set_xlabel('Years')
ax.set_ylabel('Stock [GtC]')
ax.set_title(r'$Y_0 = $' + str(y0[0]))
ax.grid()
plt.savefig('Plots/logistic_func_y0_01.png')
plt.show()

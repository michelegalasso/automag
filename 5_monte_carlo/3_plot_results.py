"""
automag.3_monte_carlo.3_plot_results.py
=======================================

Script which plots the results of the Vampire run.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit


def curve(x, Tc, beta):
    return np.where(x < Tc, np.sign(1 - x/Tc) * np.abs(1 - x/Tc) ** beta, 0)


filename = 'output'

results = np.loadtxt(filename)
temperature = results[:, 0]
magnetization = results[:, -1]

pars, cov = curve_fit(f=curve, xdata=temperature, ydata=magnetization, p0=[900, 0.34])

# curve for plot
data_curve = curve(temperature, pars[0], pars[1])

plt.rcParams.update({'font.size': 19})
plt.figure(figsize=(8, 6))
plt.plot(temperature, magnetization, label='Monte Carlo')
plt.plot(temperature, data_curve, label='analytical fit')
plt.xlabel('Temperature (K)')
plt.ylabel('Mean magnetization length')
plt.legend()
# plt.show()
plt.savefig('magnetization.png')

print(f'Fitted critical exponent beta = {pars[1]:.2f}')
print(f'Fitted critical temperature Tc = {pars[0]:.0f} K')

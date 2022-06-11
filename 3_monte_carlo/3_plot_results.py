"""
automag.3_monte_carlo.3_plot_results.py
=======================================

Script which plots the results of the Vampire run.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

import numpy as np
import matplotlib.pyplot as plt


def base(x, a):
    return 1 - x/a


exponent = 0.34
filename = 'output'

results = np.loadtxt(filename)
temperature = results[:, 0]
magnetization = results[:, -1]

tc = 1
trial = 1
previous_residual = np.inf
while trial < max(temperature):
    bloch_curve = np.array([base(item, trial) ** exponent if base(item, trial) > 0 else 0 for item in temperature])
    residual = np.sum(np.subtract(magnetization, bloch_curve) ** 2)
    if residual < previous_residual:
        tc = trial
        previous_residual = residual
    trial += 1

# Bloch curve for plot
bloch_curve = np.array([base(item, tc) ** exponent if base(item, tc) > 0 else 0 for item in temperature])

plt.rcParams.update({'font.size': 19})
plt.figure(figsize=(8, 6))
plt.plot(temperature, magnetization, label='Monte Carlo')
plt.plot(temperature, bloch_curve, label='theory')
plt.xlabel('Temperature (K)')
plt.ylabel('Mean magnetization length')
plt.legend()
# plt.show()
plt.savefig('magnetization.png')

print(f'The estimated critical temperature is {tc} K.')

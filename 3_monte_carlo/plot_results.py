"""
automag.3_monte_carlo.plot_results.py
=====================================

Script which plots the results of the Vampire run.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

import numpy as np
import matplotlib.pyplot as plt

RESULTS = 'output'

results = np.loadtxt(RESULTS)
temperature = results[:, 0]
magnetization = results[:, -1]

plt.rcParams.update({'font.size': 19})
plt.figure(figsize=(8, 6))
plt.plot(temperature, magnetization)
plt.xlabel('Temperature (K)')
plt.ylabel('Mean magnetization length')
# plt.show()
plt.savefig('magnetization.png')

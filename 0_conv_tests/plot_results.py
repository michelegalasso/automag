"""
automag.0_conv_tests.plot_results
==================================

Script which plots results of convergence tests.

.. codeauthor:: Michele Galasso <michele.galasso@skoltech.ru>
"""

import numpy as np
import matplotlib.pyplot as plt

# input parameters: just comment what you don't need
MODE = 'encut'
# MODE = 'kgrid'
ATOMS_IN_CELL = 30
XMIN = 650
XMAX = 850


def single_plot(X, Y, label=None):
    X = np.array(X, dtype=int)
    Y = np.array(Y, dtype=float) / ATOMS_IN_CELL

    # plot only between XMIN and XMAX
    if 'XMAX' in globals() and 'XMIN' in globals():
        indices = np.where(np.logical_and(X >= XMIN, X <= XMAX))
        X = X[indices]
        Y = Y[indices]

    # plot
    if label is None:
        plt.plot(X, Y, 'o')
    else:
        # sort the arrays
        indices = np.argsort(X)
        X = X[indices]
        Y = Y[indices]

        plt.plot(X, Y, 'o-', label=label)


# set figure size
plt.figure(figsize=(16, 9))

# read only lines which do not start with space
lines = []
with open(f'Fe12O18_{MODE}.txt', 'r') as f:
    for line in f:
        if line[0] != ' ':
            lines.append(line)

# extract the results
params, energies, = [], []
for line in lines:
    values = line.split()
    params.append(values[0].strip(MODE))
    energies.append(values[-1].split('=')[1])

if MODE == 'encut':
    single_plot(params, energies)
elif MODE == 'kgrid':
    sigmas, kptss = [], []
    for param in params:
        sigma, kpts = param.split('-')
        sigmas.append(sigma)
        kptss.append(kpts)

    sigmas = np.array(sigmas)
    kptss = np.array(kptss)
    energies = np.array(energies)
    sigma_values = set(sigmas)
    for sigma_value in sigma_values:
        indices = np.where(sigmas == sigma_value)
        single_plot(kptss[indices], energies[indices], label=sigma_value)
else:
    raise ValueError('MODEs currently implemented: encut, kgrid.')

plt.legend()
plt.savefig(f'{MODE}.png')
# plt.show()
"""
automag.0_conv_tests.2_plot_results
===================================

Script which plots results of convergence tests.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

from input import mode, poscar_file

import os
import numpy as np
import matplotlib.pyplot as plt

from ase.io import read

# increase matplotlib pyplot font size
plt.rcParams.update({'font.size': 20})

# create an ase atoms object
path_to_poscar = '../geometries/' + poscar_file
atoms = read(path_to_poscar)

# minimum and maximum values on the x axis for plotting, can be omitted
# XMIN = 250
# XMAX = 1000


def single_plot(X, Y, label=None):
    X = np.array(X, dtype=int)
    Y = np.array(Y, dtype=float) / len(atoms)

    # plot only between XMIN and XMAX
    if 'XMAX' in globals() and 'XMIN' in globals():
        indices = np.where(np.logical_and(X >= XMIN, X <= XMAX))
        X = X[indices]
        Y = Y[indices]

    # sort the arrays
    indices = np.argsort(X)
    X = X[indices]
    Y = Y[indices]

    for x, y in zip(X[:-1], Y[:-1]):
        if abs(y - Y[-1]) < 0.001:
            if mode == 'encut':
                print(f'ENCUT = {x} eV gives an error of less than 1 meV/atom w. r. t. the most accurate result.')
            else:
                print(f'kpts = {x} eV with {label} gives an error of less than 1 meV/atom w. r. t. '
                      f'the most accurate result.')
            break

    # plot
    if label is None:
        plt.plot(X, Y, 'o')
    else:
        plt.plot(X, Y, 'o-', label=label)


# set figure size
plt.figure(figsize=(16, 9))

# read only lines which do not start with space
lines = []
calcfold = os.path.join(os.environ.get('AUTOMAG_PATH'), 'CalcFold')
with open(os.path.join(calcfold, f"{atoms.get_chemical_formula(mode='metal', empirical=True)}_{mode}.txt"), 'r') as f:
    for line in f:
        if line[0] != ' ':
            lines.append(line)

# extract the results
params, energies, = [], []
for line in lines:
    values = line.split()
    params.append(values[0].strip(mode))
    energies.append(values[-1].split('=')[1])

if mode == 'encut':
    single_plot(params, energies)

elif mode == 'kgrid':
    sigmas, kptss = [], []
    for param in params:
        sigma, kpts = param.split('-')
        sigmas.append(sigma)
        kptss.append(kpts)

    sigmas = np.array(sigmas)
    kptss = np.array(kptss)
    energies = np.array(energies)
    sigma_values = sorted(set(sigmas))
    for sigma_value in sigma_values:
        indices = np.where(sigmas == sigma_value)
        single_plot(kptss[indices], energies[indices], label=f'SIGMA = {sigma_value}')
else:
    raise ValueError('MODEs currently implemented: encut, kgrid.')

if mode == 'encut':
    plt.xlabel('ENCUT (eV)')
else:
    plt.xlabel(r'$R_k$')
    plt.legend()

plt.ylabel('energy (eV/atom)')
plt.savefig(f'{mode}.png', bbox_inches='tight')
# plt.show()

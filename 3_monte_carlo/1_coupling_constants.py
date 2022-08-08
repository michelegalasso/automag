"""
automag.3_monte_carlo.1_coupling_constants.py
=============================================

Script which computes the coupling constants between magnetic atoms.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

from input import *

import os
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from pymatgen.core.structure import Structure

# initialize variables
structure = None
states = None
energies = None

# read input files from previous step
for item in os.listdir('../2_coll'):
    rel_path = os.path.join('../2_coll', item)
    if os.path.isfile(rel_path):
        if item.startswith('setting') and item.endswith('.vasp'):
            structure = Structure.from_file(rel_path)
        if item.startswith('states') and item.endswith('.txt'):
            with open(rel_path, 'rt') as f:
                states = json.load(f)
        if item.startswith('energies') and item.endswith('.txt'):
            with open(rel_path, 'rt') as f:
                energies = json.load(f)

if structure is None:
    raise IOError('No setting file found in ../2_coll folder.')
if states is None:
    raise IOError('No states file found in ../2_coll folder.')
if energies is None:
    raise IOError('No energies file found in ../2_coll folder.')

# find out which atoms are magnetic
for element in structure.composition.elements:
    if 'magnetic_atoms' not in globals():
        element.is_magnetic = element.is_transition_metal
    else:
        if element.name in magnetic_atoms:
            element.is_magnetic = True
        else:
            element.is_magnetic = False

# from eV/atom to total energy of the unit cell
energies = structure.num_sites * np.array(energies)


def system(configurations):
    matrix = []
    for item in configurations:
        equation = [1]
        for distance in unique_distances:
            count = 0
            for atom1, atom2, d in zip(center_indices, point_indices, distances):
                if np.isclose(d, distance, atol=0.02):
                    count += item[atom1] * item[atom2]
            equation.append(-count // 2)
        matrix.append(equation)
    return np.array(matrix)


non_magnetic_atoms = [element.symbol for element in structure.composition.elements if not element.is_magnetic]
structure.remove_species(non_magnetic_atoms)
center_indices, point_indices, offset_vectors, distances = structure.get_neighbor_list(cutoff_radius)

# get unique distances
unique_distances, counts = np.unique(np.around(distances, 2), return_counts=True)

# create fit and control group
fit_group_size = 1 - control_group_size
index = int(round(len(states) * fit_group_size))
configurations_fit = np.array(states[:index])
configurations_control = np.array(states[index:])
energies_fit = energies[:index]
energies_control = energies[index:]

A = system(configurations_fit)
values = np.linalg.lstsq(A, energies_fit, rcond=None)

# DEBUG
# all_states = [[1]]
# for _ in range(structure.num_sites - 1):
#     new = []
#     for item in all_states:
#         new.append(item + [1])
#         new.append(item + [-1])
#     all_states = new
#
# tmp = system(all_states)
# pred = tmp @ values[0]
# pred -= min(pred)
# pred *= 1000      go to meV/unit cell
# pred /= 40        go to meV/atom

B = system(configurations_control)
predictions = B @ values[0]
PCC = np.corrcoef(predictions, energies_control)

plt.rcParams.update({'font.size': 13})
plt.locator_params(axis='x', nbins=5)
plt.xlabel(f'Heisenberg model energy (eV)')
plt.ylabel(f'DFT energy (eV)')

# print results
if np.linalg.matrix_rank(A) == len(unique_distances) + 1:
    coupling_constants = values[0][1:] * 1.60218e-19
    print(f'distances between neighbors: {unique_distances.tolist()}')
    print(f'counts: {counts.tolist()}')
    print(f'coupling constants: {np.array2string(coupling_constants, precision=8, separator=", ")}')

    if append_coupling_constants:
        with open('input.py', 'a') as f:
            f.write('\n# LINE ADDED BY THE SCRIPT 1_coupling_constants.py')
            f.write(f'\ndistances_between_neighbors = {unique_distances.tolist()}\n')

            f.write('\n# LINE ADDED BY THE SCRIPT 1_coupling_constants.py')
            f.write(f'\ncoupling_constants = {np.array2string(coupling_constants, precision=8, separator=", ")}\n')

    plt.scatter(predictions, energies_control, label=f'PCC: {PCC[0, 1]:.2f}')
    ax = plt.gca()
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
    plt.legend()
    # plt.show()
    plt.savefig('model.png', bbox_inches='tight')
    print(f'PCC: {PCC[0, 1]:.2f}')
else:
    print(f'ERROR: SYSTEM OF {np.linalg.matrix_rank(A)} INDEPENDENT EQUATION(S) '
          f'IN {len(unique_distances) + 1} UNKNOWNS!')

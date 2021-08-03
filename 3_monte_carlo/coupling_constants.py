"""
automag.3_monte_carlo.coupling_constants.py
===========================================

Script which computes the coupling constants between magnetic atoms.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

import json
import numpy as np
import matplotlib.pyplot as plt

from pymatgen.core.structure import Structure

# input parameters
GEOMETRY = '../input/Fe2O3-alpha.vasp'
NON_MAGNETIC_ATOMS = ['O']
CUTOFF_RADIUS = 3.8
N_MAGNETIC_SITES = 12


def system(configurations):
    matrix = []
    for item in configurations:
        equation = [1]
        for distance in unique_distances:
            count = 0
            for atom1, atom2, d in zip(center_indices, point_indices, distances):
                if np.isclose(d, distance, atol=0.05):
                    count += item[atom1] * item[atom2]
            equation.append(-count // 2)
        matrix.append(equation)
    return np.array(matrix)


structure = Structure.from_file(GEOMETRY)
structure.remove_species(NON_MAGNETIC_ATOMS)
center_indices, point_indices, offset_vectors, distances = structure.get_neighbor_list(CUTOFF_RADIUS)

# get unique distances
unique_distances, counts = np.unique(np.around(distances, 2), return_counts=True)

# read configurations from file
with open('configurations.txt', 'rt') as f:
    configurations = json.load(f)

# read enthalpies from file
with open('energies.txt', 'rt') as f:
    energies = json.load(f)

# create fit and control group 60% to 40%
index = int(round(len(configurations) * 0.6))
configurations_fit = np.array(configurations[:index])
configurations_control = np.array(configurations[index:])
energies_fit = np.array(energies[:index])
energies_control = np.array(energies[index:])

A = system(configurations_fit)
values = np.linalg.lstsq(A, energies_fit, rcond=None)

B = system(configurations_control)
predictions = B @ values[0]
PCC = np.corrcoef(predictions, energies_control)

magnetic_atom = structure
plt.rcParams.update({'font.size': 19})
plt.figure(figsize=(16, 9))
plt.xlabel(f'Heisemberg model free energy (eV/{structure.types_of_species[0].name} atom)')
plt.ylabel(f'DFT free energy (eV/{structure.types_of_species[0].name} atom)')

# print results
if np.linalg.matrix_rank(A) == len(unique_distances) + 1:
    print(f'distances between neighbours: {unique_distances.tolist()}')
    print(f'counts: {counts.tolist()}')
    print(f'coupling constants: {values[0][1:] * 1.60218e-19}')
    plt.scatter(predictions / N_MAGNETIC_SITES, energies_control / N_MAGNETIC_SITES, label=f'PCC: {PCC[0, 1]:.2f}')
    plt.legend()
    # plt.show()
    plt.savefig('model.png')
else:
    print(f'ERROR: SYSTEM OF {np.linalg.matrix_rank(A)} INDEPENDENT EQUATION(S) '
          f'IN {len(unique_distances) + 1} UNKNOWNS!')

"""
automag.3_monte_carlo.coupling_constants.py
===========================================

Script which computes the coupling constants between magnetic atoms.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

import json
import numpy as np

from pymatgen.core.structure import Structure

# input parameters
RELAXED_GEOMETRY = 'Fe2O3-alpha-afm1-relaxed.vasp'
NON_MAGNETIC_ATOMS = ['O']
NEIGHBORS_RADIUS = 4.0

structure = Structure.from_file(RELAXED_GEOMETRY)
structure.remove_species(NON_MAGNETIC_ATOMS)
center_indices, point_indices, offset_vectors, distances = structure.get_neighbor_list(NEIGHBORS_RADIUS)

# get unique distances
unique_distances, counts = np.unique(np.around(distances, 2), return_counts=True)
print(f'distances between neighbours: {unique_distances.tolist()}')
print(f'counts: {counts.tolist()}')

# read configurations from file
with open('configurations.txt', 'rt') as f:
    configurations = json.load(f)

# read enthalpies from file
with open('enthalpies.txt', 'rt') as f:
    b = json.load(f)

A = []
for item in configurations:
    equation = [1]
    for distance in unique_distances:
        count = 0
        for atom1, atom2, d in zip(center_indices, point_indices, distances):
            if np.isclose(d, distance, atol=0.05):
                count += item[atom1] * item[atom2]
        equation.append(-count // 2)
    A.append(equation)
A = np.array(A)

values = np.linalg.lstsq(A, b, rcond=None)
print(f'coupling constants: {values[0][1:] * 1.60218e-19}')

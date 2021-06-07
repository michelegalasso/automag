"""
automag.3_monte_carlo.run_vampire.py
====================================

Script which writes the unit cell file for running Vampire.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""
import numpy as np
from pymatgen.core.structure import Structure

# input parameters
RELAXED_GEOMETRY = 'La2Nd2H40-fm-relaxed.vasp'
NON_MAGNETIC_ATOMS = ['La', 'H']
CUTOFF_RADIUS = 4.0
DISTANCES_BETWEEN_NEIGHBORS = [3.57]
COUPLING_CONSTANTS = [4.25219829e-18]
NUM_MATERIALS = 1           # set 1 for Curie temperature, 2 for Néel temperature

structure = Structure.from_file(RELAXED_GEOMETRY)
structure.remove_species(NON_MAGNETIC_ATOMS)
center_indices, point_indices, offset_vectors, distances = structure.get_neighbor_list(CUTOFF_RADIUS)

# get thresholds for comparing atomic distances
thresholds = [(a + b) / 2 for a, b in zip(DISTANCES_BETWEEN_NEIGHBORS[:-1], DISTANCES_BETWEEN_NEIGHBORS[1:])]
thresholds = [0.0] + thresholds + [100.0]

# print unit cell size for VAMPIRE input
with open('vamp.ucf', 'w') as f:
    f.write('# Unit cell size:\n')
    f.write(f'{structure.lattice.a:19.16f}  {structure.lattice.b:19.16f}  {structure.lattice.c:19.16f}\n')

    f.write('# Unit cell vectors:\n')
    for i, (vector, modulus) in enumerate(zip(structure.lattice.matrix, structure.lattice.abc)):
        v = vector / modulus
        f.write(f'{v[0]: 10.8e}   {v[1]: 10.8e}   {v[2]: 10.8e}\n')

    # print fractional coordinates for VAMPIRE input
    f.write('# Atoms num_atoms num_materials; id cx cy cz mat cat hcat\n')
    f.write(f'{len(structure):d} {NUM_MATERIALS}\n')
    for i, coord in enumerate(structure.frac_coords):
        for j, item in enumerate(coord):
            if item < 0:
                coord[j] += 1
            elif item >= 1:
                coord[j] -= 1
        material = i // (len(structure) // NUM_MATERIALS)
        f.write(f'{i:2d}   {coord[0]:18.16f}  {coord[1]:18.16f}  {coord[2]:18.16f}  {material} 0 0\n')

    f.write('# Interactions n exctype; id i j dx dy dz Jij\n')
    f.write(f'{len(center_indices)} isotropic\n')
    for i, (atom1, atom2, offset, distance) in enumerate(zip(center_indices, point_indices, offset_vectors, distances)):
        for low_lim, high_lim, coupling_constant in zip(thresholds[:-1], thresholds[1:], COUPLING_CONSTANTS):
            if low_lim < distance < high_lim:
                f.write(f'{i:3d}   {atom1:2d}  {atom2:2d}  '
                        f'{int(offset[0]):2d} {int(offset[1]):2d} {int(offset[2]):2d}   {coupling_constant: 6.4e}\n')

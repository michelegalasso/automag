"""
automag.3_monte_carlo.2_run_vampire.py
======================================

Script which runs Vampire for Monte Carlo simulation.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

from input import *

import os
import numpy as np

from pymatgen.core.structure import Structure

# raise IOError if no trials folder
if not os.path.isdir('../2_coll/trials'):
    raise IOError('No trials folder found in ../2_coll.')

setting = 1
magmom = None
while os.path.isfile(f'../2_coll/trials/configurations{setting:03d}.txt'):
    with open(f'../2_coll/trials/configurations{setting:03d}.txt', 'rt') as f:
        for line in f:
            values = line.split()
            if values[0] == configuration:
                magmom = [int(item) for item in values[1:]]
                break
        if magmom is not None:
            break
    setting += 1

if magmom is None:
    raise IOError(f'configuration {configuration} not found.')

# get materials from magmom
materials = [0 if item >= 0 else 1 for item in magmom]

# create a pymatgen Structure object
structure = Structure.from_file(f'../2_coll/trials/setting{setting:03d}.vasp')

# find out which atoms are magnetic
for element in structure.composition.elements:
    if 'magnetic_atoms' not in globals():
        element.is_magnetic = element.is_transition_metal
    else:
        if element.name in magnetic_atoms:
            element.is_magnetic = True
        else:
            element.is_magnetic = False

non_magnetic_atoms = [element.symbol for element in structure.composition.elements if not element.is_magnetic]
structure.remove_species(non_magnetic_atoms)
center_indices, point_indices, offset_vectors, distances = structure.get_neighbor_list(cutoff_radius)

# get thresholds for comparing atomic distances
thresholds = [(a + b) / 2 for a, b in zip(distances_between_neighbors[:-1], distances_between_neighbors[1:])]
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
    f.write(f'{len(structure):d} {len(np.unique(materials))}\n')
    for i, (coord, material) in enumerate(zip(structure.frac_coords, materials)):
        coord += 0.     # gets rid of the minus zero values
        for j, item in enumerate(coord):
            if item < 0:
                coord[j] += 1
            elif item >= 1:
                coord[j] -= 1
        f.write(f'{i:2d}   {coord[0]:18.16f}  {coord[1]:18.16f}  {coord[2]:18.16f}  {material} 0 0\n')

    f.write('# Interactions n exctype; id i j dx dy dz Jij\n')
    f.write(f'{len(center_indices)} isotropic\n')
    for i, (atom1, atom2, offset, distance) in enumerate(zip(center_indices, point_indices, offset_vectors, distances)):
        for low_lim, high_lim, coupling_constant in zip(thresholds[:-1], thresholds[1:], coupling_constants):
            if low_lim < distance < high_lim:
                f.write(f'{i:3d}   {atom1:2d}  {atom2:2d}  '
                        f'{int(offset[0]):2d} {int(offset[1]):2d} {int(offset[2]):2d}   {coupling_constant: 6.4e}\n')

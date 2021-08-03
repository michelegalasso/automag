"""
automag.3_monte_carlo.get_configurations.py
===========================================

Script which computes the magnetic configurations associates with different poscar files.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

import os
import json
import numpy as np

from pymatgen.core.structure import Structure
from pymatgen.core.composition import Composition

from common.insert_elements_in_poscar import insert_elements_in_poscar

# input parameters
COMPOSITION = Composition('Fe12O18')
MAGNETIC_ATOM = 'Fe'
ENUMLIB_DIR = '../2_coll/enumlib'
PREVIOUS_RESULTS = '../2_coll/Fe12O18_singlepoint.txt'

# reference magnetic configurations
fm_reference = np.array(12 * [4.0] + 18 * [0.0])                    # FM configuration
afm_reference = np.array(6 * [4.0] + 6 * [-4.0] + 18 * [0.0])       # AFM configuration


def get_final_magmoms_and_energy(configuration, previous_results):
    ind_1 = previous_results.index(configuration + ' ')
    ind_2 = previous_results.index('\n', ind_1)
    energy = float(previous_results[ind_1:ind_2].split()[-1].split('=')[1])

    ind_3 = previous_results.index('final_magmoms=[', ind_2) + 15
    ind_4 = previous_results.index(']', ind_3)
    final_magmoms = [float(item) for item in previous_results[ind_3:ind_4].split()]

    return final_magmoms, energy


# select filenames starting with the string 'vasp'
filenames = []
for filename in os.listdir(ENUMLIB_DIR):
    if filename.startswith('vasp'):
        filenames.append(filename)
filenames.sort(key=lambda x: int(x.split('.')[1]))

# read previous results
with open(PREVIOUS_RESULTS, 'rt') as f:
    previous_results = f.read()

# work out configurations
configurations = []
energies = []
geometries = []
previous_poscar = None
for filename in filenames:
    # read current poscar
    with open(os.path.join(ENUMLIB_DIR, filename), 'r') as f:
        current_poscar = f.read()

    config_number = filename.split('.')[1]
    current_poscar = insert_elements_in_poscar(current_poscar, COMPOSITION, MAGNETIC_ATOM)

    if 'afm' + config_number in previous_results:
        final_magmoms, energy = get_final_magmoms_and_energy('afm' + config_number, previous_results)
    else:
        continue

    if previous_poscar:
        if previous_poscar.split('\n')[2:5] == current_poscar.split('\n')[2:5]:
            geometries[-1].append(current_poscar)

            if np.array_equal(np.around(final_magmoms), afm_reference) or \
                    np.array_equal(np.around(final_magmoms), -afm_reference):
                current_structure = Structure.from_str(current_poscar, fmt='poscar')

                mag_indices = current_structure.indices_from_symbol(MAGNETIC_ATOM)
                current_coords = np.array(current_poscar.split('\n'))[np.array(mag_indices) + 8]
                reference_coords = np.array(geometries[-1][0].split('\n'))[np.array(mag_indices) + 8]

                current_coords = current_coords.tolist()
                reference_coords = reference_coords.tolist()
                indices = [current_coords.index(coord) for coord in reference_coords]

                configuration = np.sign(afm_reference).astype(int)
                configurations[-1].append(configuration[indices])
                energies[-1].append(energy)
        else:
            previous_poscar = None
            geometries.append([current_poscar])

            if np.array_equal(np.around(final_magmoms), afm_reference) or \
                    np.array_equal(np.around(final_magmoms), -afm_reference):
                configuration = np.sign(afm_reference).astype(int)
                configurations.append([configuration[configuration != 0]])
                energies.append([energy])
    else:
        geometries.append([current_poscar])

        if np.array_equal(np.around(final_magmoms), afm_reference) or \
                np.array_equal(np.around(final_magmoms), -afm_reference):
            configuration = np.sign(afm_reference).astype(int)
            configurations.append([configuration[configuration != 0]])
            energies.append([energy])

        if config_number == '1':
            final_magmoms, energy = get_final_magmoms_and_energy('fm', previous_results)
            if np.array_equal(np.around(final_magmoms), fm_reference) or \
                    np.array_equal(np.around(final_magmoms), -fm_reference):
                configuration = np.sign(fm_reference).astype(int)
                configurations[-1].append(configuration[configuration != 0])
                energies[-1].append(energy)

    previous_poscar = current_poscar

# select geometry with most configurations
geometry_index = 0
for i, l in enumerate(configurations):
    if len(l) > len(configurations[geometry_index]):
        geometry_index = i

# print geometry with most configurations and its abundance
print('Geometry: afm{}'.format(geometries[geometry_index][0].split('\n')[0].split()[-1]))
print('Abundance: {}'.format(len(configurations[geometry_index])))

# write configurations to file
with open('configurations.txt', 'wt') as f:
    json.dump([a.tolist() for a in configurations[geometry_index]], f)

# write energies to file
with open('energies.txt', 'wt') as f:
    json.dump(energies[geometry_index], f)

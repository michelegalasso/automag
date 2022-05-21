"""
automag.2_coll.1_submit
=======================

Script which runs enumlib and submits calculations.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

from input import *

import os
import subprocess
import numpy as np

from itertools import product
from pymatgen.io.vasp import Poscar
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from common.SubmitFirework import SubmitFirework


def launch_enumlib(count, split):
    os.mkdir(f'enumlib{count}')
    os.chdir(f'enumlib{count}')

    with open('struct_enum.in', 'w') as f:
        f.write('generated by Automag\n')
        f.write('bulk\n')

        for lat_vector in symmetrized_structure.lattice.matrix:
            for component in lat_vector:
                f.write(f'{component:14.10f}        ')
            f.write('\n')

        case = len(split) + sum(split)
        f.write(f'  {case} -nary case\n')
        f.write(f'    {symmetrized_structure.num_sites} # Number of points in the multilattice\n')

        offset = 0
        for i, (s, wyckoff) in enumerate(zip(split, symmetrized_structure.equivalent_sites)):
            if s:
                offset += 1

            for atom in wyckoff:
                for component in atom.coords:
                    f.write(f'{component:14.10f}        ')
                if s:
                    f.write(f'{i + offset - 1}/{i + offset}\n')
                else:
                    f.write(f'{i + offset}\n')

        f.write(f'    1 {supercell_size}   # Starting and ending cell sizes for search\n')
        f.write('0.10000000E-06 # Epsilon (finite precision parameter)\n')
        f.write('full list of labelings\n')
        f.write('# Concentration restrictions\n')

        for s, wyckoff in zip(split, symmetrized_structure.equivalent_sites):
            if s:
                for _ in range(2):
                    f.write(f'{len(wyckoff):4d}')
                    f.write(f'{len(wyckoff):4d}')
                    f.write(f'{symmetrized_structure.num_sites * 2:4d}\n')
            else:
                f.write(f'{len(wyckoff) * 2:4d}')
                f.write(f'{len(wyckoff) * 2:4d}')
                f.write(f'{symmetrized_structure.num_sites * 2:4d}\n')

    # process = subprocess.Popen('/home/michele/softs/enumlib/src/enum.x')
    process = subprocess.Popen('enum.x')
    try:
        process.wait(timeout=60)
    except subprocess.TimeoutExpired:
        process.kill()

    # os.system('/home/michele/softs/enumlib/aux_src/makeStr.py 1 500')
    os.system('makeStr.py 1 500')

    for j in range(501):
        if os.path.isfile(f'vasp.{j + 1}'):
            conf_poscar = Poscar.from_file(f'vasp.{j + 1}')
        else:
            break

        if conf_poscar.structure.lattice in lattices:
            index = lattices.index(conf_poscar.structure.lattice)
            current_coords = conf_poscar.structure.frac_coords.tolist()
            reference_coords = coordinates[index].tolist()
            mapping = [current_coords.index(coord) for coord in reference_coords]
        else:
            lattices.append(conf_poscar.structure.lattice)
            coordinates.append(conf_poscar.structure.frac_coords)
            configurations.append([])
            index = len(lattices) - 1
            mapping = list(range(len(conf_poscar.structure.frac_coords)))

        k = 0
        groups = []
        site_magmoms = []
        for s, states in zip(split, wyckoff_magmoms):
            if s:
                groups.append([k, k + 1])
                for _ in range(2):
                    site_magmoms.append([spin_value, -spin_value])
                    k += 1
            else:
                site_magmoms.append(states)
                k += 1

        # use a flag to check that all sites which have been split are AFM
        for conf in product(*site_magmoms):
            flag = True
            conf_array = np.array(conf)
            for group in groups:
                if sum(conf_array[group]) != 0:
                    flag = False
            if flag:
                configuration = np.repeat(conf, conf_poscar.natoms)
                transformed_configuration = configuration[mapping]
                if transformed_configuration.tolist() not in configurations[index] and \
                        [-item for item in transformed_configuration] not in configurations[index] and \
                        sum(abs(transformed_configuration)) != 0:
                    configurations[index].append(transformed_configuration.tolist())

    os.chdir('..')


# full path to poscar file
path_to_poscar = '../geometries/' + poscar_file

# create Structure and SymmetrizedStructure objects
structure = Structure.from_file(path_to_poscar)
analyzer = SpacegroupAnalyzer(structure)
symmetrized_structure = analyzer.get_symmetrized_structure()

# find out which atoms are magnetic
for element in structure.composition.elements:
    if 'magnetic_atoms' not in globals():
        element.is_magnetic = element.is_transition_metal
    else:
        if element.name in magnetic_atoms:
            element.is_magnetic = True
        else:
            element.is_magnetic = False

if os.path.exists('trials'):
    print('Cannot create a folder named trials: an object with the same name already exists.')
    exit()

os.mkdir('trials')
os.chdir('trials')

# geometrical settings and respective lists of magnetic configurations
lattices = [structure.lattice]
coordinates = [structure.frac_coords]
configurations = [[]]

# get the multiplicities of each Wyckoff position
multiplicities = [len(item) for item in symmetrized_structure.equivalent_indices]

wyckoff_magmoms = []
for i, wyckoff in enumerate(symmetrized_structure.equivalent_sites):
    if wyckoff[0].specie.is_magnetic:
        wyckoff_magmoms.append([spin_value, 0, -spin_value])
    else:
        wyckoff_magmoms.append([0])

# get all possible configurations without any splitting
for conf in product(*wyckoff_magmoms):
    configuration = np.repeat(conf, multiplicities)
    if [-item for item in configuration] not in configurations[0] and sum(abs(configuration)) != 0:
        configurations[0].append(configuration.tolist())

# split all possible combinations of Wyckoff positions
splits = []
possibilities = [[0, 1] if len(item) > 1 else [0] for item in wyckoff_magmoms]
for split in product(*possibilities):
    if sum(split) != 0:
        splits.append(split)

for i, split in enumerate(splits):
    launch_enumlib(i + 1, split)

# merge the first and second settings if they are equal
if len(coordinates) > 1 and len(coordinates[0]) == len(coordinates[1]):
    del lattices[0]
    del coordinates[0]
    configurations[0].extend(configurations[1])
    del configurations[1]

# write output and submit calculations
afm_count = 1
fim_count = 1
original_ch_symbols = [atom.name for atom in structure.species]
for i, (lattice, frac_coords, confs) in enumerate(zip(lattices, coordinates, configurations)):
    magnification = len(frac_coords) // len(structure.frac_coords)
    ch_symbols = np.repeat(original_ch_symbols, magnification)
    setting = Structure(lattice, ch_symbols, frac_coords)
    setting.to(fmt='poscar', filename=f'setting{i + 1:03d}.vasp')
    mask = [item.is_magnetic for item in setting.species]

    for conf in confs:
        conf_array = np.array(conf)
        with open(f'configurations{i + 1:03d}.txt', 'a') as f:
            if np.sum(conf) == spin_value * np.sum(mask):
                state = 'fm'
            elif np.sum(np.abs(conf)) == 0:
                state = 'nm'
            elif np.sum(conf) == 0:
                state = 'afm' + str(afm_count)
                afm_count += 1
            else:
                state = 'fim' + str(fim_count)
                fim_count += 1

            f.write(f'{state:>6s}  ')
            f.write(' '.join(f'{e:2d}' for e in conf_array[mask]))
            f.write('\n')

        run = SubmitFirework(f'setting{i + 1:03d}.vasp', mode='singlepoint', fix_params=params, magmoms=conf,
                             name=state)
        run.submit()

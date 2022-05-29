"""
automag.2_coll.2_plot_results
=============================

Script which plots results of magnetic relaxations.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

from input import *

import os
import json
import shutil
import numpy as np
import matplotlib.pyplot as plt

from copy import copy
from ase.io import read
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# take care of the case when lower_cutoff has not been specified
if 'lower_cutoff' not in globals():
    lower_cutoff = 0

# increase matplotlib pyplot font size
plt.rcParams.update({'font.size': 20})

# set figure size
plt.figure(figsize=(16, 9))

# create an ase atoms object
path_to_poscar = '../geometries/' + poscar_file
atoms = read(path_to_poscar)

# get the multiplicities of each Wyckoff position
structure = Structure.from_file(path_to_poscar)
analyzer = SpacegroupAnalyzer(structure)
symmetrized_structure = analyzer.get_symmetrized_structure()
multiplicities = np.array([len(item) for item in symmetrized_structure.equivalent_indices])

# path to results file with data
calcfold = os.path.join(os.environ.get('AUTOMAG_PATH'), 'CalcFold')
output_file = os.path.join(calcfold, f"{atoms.get_chemical_formula(mode='metal', empirical=True)}_singlepoint.txt")

# exit if no trials folder
if not os.path.isdir('trials'):
    raise IOError('No trials folder found.')

# read lines
lines, maginfos, presents = [], [], []
with open(output_file, 'rt') as f:
    for line in f:
        if line[0] != ' ':
            presents.append(line.split()[0])
            lines.append(line)
        else:
            maginfos.append(line)

data = {}
setting = 1
while os.path.isfile(f'trials/configurations{setting:03d}.txt'):
    with open(f'trials/configurations{setting:03d}.txt', 'rt') as f:
        for line in f:
            values = line.split()
            init_state = values[0]
            if init_state in presents:
                dct = {
                    'setting': setting,
                    'init_spins': [int(item) for item in values[1:]],
                }
                data[init_state] = dct
    setting += 1

# keep the maximum value of setting
max_setting = copy(setting)

# extract the results
all_final_magmoms = []
for line, maginfo in zip(lines, maginfos):
    values = line.split()
    init_state = values[0]

    if values[2].split('=')[1] == 'NONCONVERGED':
        print(f'WARNING: energy calculation of {values[0]} did not converge and will be excluded from the graph.')
        del data[init_state]
    else:
        initial, final = maginfo.split('final_magmoms=')
        initial = initial[initial.index('[') + 1:initial.index(']')]
        final = final[final.index('[') + 1:final.index(']')]
        initial = np.array(initial.split(), dtype=float)
        final = np.array(final.split(), dtype=float)

        all_final_magmoms.extend(final.tolist())

        # exclude low-spin configurations
        flag = False
        mask = np.nonzero(initial)
        if np.all(np.abs(final[mask]) > lower_cutoff) or values[0] == 'nm':
            flag = True

        magnification = len(initial) // sum(multiplicities)
        current_multiplicities = magnification * multiplicities

        start = 0
        kept = True
        for multiplicity in multiplicities:
            w_initial = initial[start:start + multiplicity]
            w_final = final[start:start + multiplicity]
            start += multiplicity

            if not np.any(w_initial):
                if np.any(np.around(w_final)):
                    kept = False
            else:
                prod = np.multiply(np.sign(w_initial), w_final)
                if prod.max() - prod.min() > 0.1:
                    kept = False

        if kept and flag:
            data[init_state]['kept_magmoms'] = True
        else:
            data[init_state]['kept_magmoms'] = False
            print(f'WARNING: {values[0]} did not keep the original magmoms and will be marked in red on the graph.')

        data[init_state]['energy'] = float(values[-1].split('=')[1]) / len(atoms)

# plt.hist(np.abs(all_final_magmoms), bins=40)
# plt.savefig('spin_distribution.png')

setting = 1
final_states = []
final_setting = setting
final_energies = []
while setting < max_setting:
    current_states = []
    current_energies = []
    for init_state, value in data.items():
        if init_state != 'nm' and value['setting'] == setting and value['kept_magmoms']:
            current_states.append(np.sign(value['init_spins']).tolist())
            current_energies.append(value['energy'])
    if len(current_states) > len(final_states):
        final_setting = setting
        final_states = current_states
        final_energies = current_energies
    setting += 1

# write states to file
with open(f'states{final_setting:03d}.txt', 'wt') as f:
    json.dump(final_states, f)

# write energies to file
with open(f'energies{final_setting:03d}.txt', 'wt') as f:
    json.dump(final_energies, f)

# copy setting file with geometry
shutil.copy(f'trials/setting{final_setting:03d}.vasp', '.')

# extract values for plot
bar_labels = []
energies = []
kept_magmoms = []
for key, value in sorted(data.items(), key=lambda t: t[0] if t[0][:2] != 'fm' else 'aaa' + t[0][2:]):
    bar_labels.append(key)
    energies.append(value['energy'])
    kept_magmoms.append(value['kept_magmoms'])

energies = np.array(energies)
kept_magmoms = np.array(kept_magmoms)

# energies from eV/magnetic atom to meV/magnetic atom
energies *= 1000

# normalize energies for plot
energies -= min(energies)
toplim = max(energies) * 1.12
energies += 0.1 * max(energies)

# plot results
repr_configurations = np.arange(1, len(energies) + 1)
plt.bar(repr_configurations[~kept_magmoms], energies[~kept_magmoms], bottom=-0.1 * max(energies), color='r')
plt.bar(repr_configurations[kept_magmoms], energies[kept_magmoms], bottom=-0.1 * max(energies), color='b')

# label bars
ax = plt.gca()
rects = ax.patches
rects.sort(key=lambda x: x.get_x())
for bar_label, rect in zip(bar_labels, rects):
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width() / 2, height - 0.09 * max(energies), bar_label,
            fontsize='small', ha='center', va='bottom', rotation='vertical')

# label axes
plt.xlabel('configurations')
plt.ylabel(f'relative energy (meV/atom)')
plt.xticks(np.arange(1, len(rects) + 1))
plt.ylim(top=toplim)

# save or show bar chart
plt.savefig('stability.png', bbox_inches='tight')
# plt.show()

print(f'The most stable configuration is {bar_labels[np.argmin(energies)]}.')
if np.logical_not(kept_magmoms[np.argmin(energies)]):
    print('WARNING: values of initial and final magnetic moments of the most stable configuration '
          'significantly differ.')

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

# take care of the case when hs_cutoff has not been specified
if 'hs_cutoff' not in globals():
    hs_cutoff = 0

# increase matplotlib pyplot font size
plt.rcParams.update({'font.size': 16})

# set figure size
plt.figure(figsize=(16, 9))

# create an ase atoms object
path_to_poscar = '../geometries/' + poscar_file
atoms = read(path_to_poscar)

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
    data[init_state]['energy'] = float(values[-1].split('=')[1]) / len(data[init_state]['init_spins'])

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
        if np.all(np.abs(final[mask]) > hs_cutoff) or values[0] == 'nm':
            flag = True

        for i, item in enumerate(final):
            if item > 0:
                final[i] = np.floor(item)
            else:
                final[i] = np.ceil(item)

        if (np.array_equal(np.sign(initial), np.sign(final)) or np.array_equal(np.sign(initial), -np.sign(final))) and \
                flag:
            data[init_state]['kept_magmoms'] = True
        else:
            data[init_state]['kept_magmoms'] = False
            print(f'WARNING: {values[0]} did not keep the original magmoms and will be marked in red on the graph.')

        data[init_state]['energy'] = float(values[-1].split('=')[1]) / len(data[init_state]['init_spins'])

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
bar_labels = [key for key in data.keys()]
energies = np.array([value['energy'] for value in data.values()])
kept_magmoms = np.array([value['kept_magmoms'] for value in data.values()])

# energies from eV/magnetic atom to meV/magnetic atom
energies *= 1000

# normalize energies for plot
energies -= min(energies)
energies += 0.1 * max(energies)

# plot results
repr_configurations = np.array(range(len(energies)))
plt.bar(repr_configurations[~kept_magmoms], energies[~kept_magmoms], bottom=-0.1 * max(energies), color='r')
plt.bar(repr_configurations[kept_magmoms], energies[kept_magmoms], bottom=-0.1 * max(energies), color='b')

# label bars
ax = plt.gca()
rects = ax.patches
rects.sort(key=lambda x: x.get_x())
for bar_label, rect in zip(bar_labels, rects):
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width() / 2, height - 0.09 * max(energies), bar_label,
            fontsize='xx-small', ha='center', va='bottom', rotation='vertical')

# label axes
plt.xlabel('configurations')
plt.ylabel(f'free energy TOTEN (meV/magnetic atom)')

# save or show bar chart
plt.savefig('stability.png')
# plt.show()

print(f'The most stable configuration is {bar_labels[np.argmin(energies)]}.')
if np.logical_not(kept_magmoms[np.argmin(energies)]):
    print('WARNING: values of initial and final magnetic moments of the most stable configuration '
          'significantly differ.')

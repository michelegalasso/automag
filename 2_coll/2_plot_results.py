"""
automag.2_coll.2_plot_results
=============================

Script which plots results of magnetic relaxations.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

from input import poscar_file

import os
import numpy as np
import matplotlib.pyplot as plt

from ase.io import read


# increase matplotlib pyplot font size
plt.rcParams.update({'font.size': 16})

# set figure size
plt.figure(figsize=(16, 9))

# create an ase atoms object
path_to_poscar = '../geometries/' + poscar_file
atoms = read(path_to_poscar)

# path to results file with data
calcfold = '../CalcFold'
data = os.path.join(calcfold, f"{atoms.get_chemical_formula(mode='metal', empirical=True)}_singlepoint.txt")

count = 0
n_magatoms = []
while os.path.isfile(f'trials/configurations{count + 1:03d}.txt'):
    with open(f'trials/configurations{count + 1:03d}.txt', 'rt') as f:
        line = f.readline()
        spins = [int(item) for item in line.split()[1:]]
    n_magatoms.append(len(spins))
    count += 1

# read lines
lines, maginfos = [], []
with open(data, 'rt') as f:
    for line in f:
        if line[0] != ' ':
            lines.append(line)
        else:
            maginfos.append(line)

    # extract the results
fw_IDs, energies, kept_magmoms = [], [], []
for line, maginfo in zip(lines, maginfos):
    values = line.split()
    n_setting = int(values[0].strip('setting'))
    fw_IDs.append(int(values[1]))
    energies.append(float(values[-1].split('=')[1]) / n_magatoms[n_setting - 1] * 1000)

    if values[2].split('=')[1] == 'NONCONVERGED':
        print(f'WARNING: energy calculation of firework {values[1]} did not converge.')

    initial, final = maginfo.split('final_magmoms=')
    initial = initial[initial.index('[') + 1:initial.index(']')]
    final = final[final.index('[') + 1:final.index(']')]
    initial = np.array(initial.split(), dtype=float)
    final = np.array(final.split(), dtype=float)

    for i, item in enumerate(final):
        if item > 0:
            final[i] = np.floor(item)
        else:
            final[i] = np.ceil(item)

    if np.array_equal(np.sign(initial), np.sign(final)) or np.array_equal(np.sign(initial), -np.sign(final)):
        kept_magmoms.append(True)
    else:
        kept_magmoms.append(False)

energies = np.array(energies)
kept_magmoms = np.array(kept_magmoms)

# sort configurations
indices = np.argsort(fw_IDs)
energies = energies[indices]
kept_magmoms = kept_magmoms[indices]

# normalize energies for plot
energies -= min(energies)
energies += 0.1 * max(energies)

# plot results
repr_configurations = np.array(range(-1, len(energies) - 1))
plt.bar(repr_configurations[~kept_magmoms], energies[~kept_magmoms], bottom=-0.1 * max(energies), color='r')
plt.bar(repr_configurations[kept_magmoms], energies[kept_magmoms], bottom=-0.1 * max(energies), color='b')

# get labels for bars
bar_labels = ['fm', 'nm']
if len(energies) > 2:
    count = 1
    for _ in energies[2:]:
        bar_labels.append(str(count))
        count += 1

# label bars
ax = plt.gca()
rects = ax.patches
rects.sort(key=lambda x: x.get_x())
for bar_label, rect in zip(bar_labels, rects):
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width() / 2, height - 0.1 * max(energies), bar_label,
            fontsize='xx-small', ha='center', va='bottom')

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

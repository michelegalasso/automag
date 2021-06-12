"""
automag.2_coll.plot_results
===========================

Script which plots results of magnetic relaxations.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

import numpy as np
import matplotlib.pyplot as plt

from pymatgen.core.composition import Composition

# input parameters
FILENAMES = ['Fe12O18_singlepoint.txt']
MAGNETIC_ATOM = 'Fe'

# increase matplotlib pyplot font size
plt.rcParams.update({'font.size': 16})

# set figure size
plt.figure(figsize=(16, 9))


def sort_configurations(conf: str):
    if conf.startswith('nm'):
        return -1
    if conf.startswith('fm'):
        return 0
    if conf.startswith('afm'):
        return int(conf.strip('afm'))


# read lines
for filename in FILENAMES:
    lines, maginfos = [], []
    with open(filename, 'r') as f:
        for line in f:
            if line[0] != ' ':
                lines.append(line)
            else:
                maginfos.append(line)

    # extract the results
    composition = Composition(filename.split('_')[0])
    configurations, energies, kept_magmoms = [], [], []
    for line, maginfo in zip(lines, maginfos):
        values = line.split()
        configurations.append(values[0])
        energies.append(float(values[-1].split('=')[1]) / composition[MAGNETIC_ATOM] * 1000)

        if values[2].split('=')[1] == 'NONCONVERGED':
            print(f'WARNING: energy calculation of {values[0]} did not converge.')

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

configurations = np.array(configurations)
energies = np.array(energies)
kept_magmoms = np.array(kept_magmoms)

# sort configurations
repr_configurations = np.array([sort_configurations(c) for c in configurations])
indices = np.argsort(repr_configurations)
repr_configurations = repr_configurations[indices]
configurations = configurations[indices]
energies = energies[indices]
kept_magmoms = kept_magmoms[indices]

# normalize energies for plot
energies -= min(energies)
energies += 0.1 * max(energies)

# plot results
plt.bar(repr_configurations[~kept_magmoms], energies[~kept_magmoms], bottom=-0.1 * max(energies), color='r')
plt.bar(repr_configurations[kept_magmoms], energies[kept_magmoms], bottom=-0.1 * max(energies), color='b')

# label bars
ax = plt.gca()
rects = ax.patches
rects.sort(key=lambda x: x.get_x())
for i, rect in enumerate(rects):
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width() / 2, height - 0.1 * max(energies), str(i),
            fontsize='xx-small', ha='center', va='bottom')

# label axes
plt.xlabel('configurations')
plt.ylabel(f'free energy TOTEN (meV/{MAGNETIC_ATOM} atom)')

# save or show bar chart
plt.savefig('stability.png')
# plt.show()

print(f'The most stable configuration is {configurations[np.argmin(energies)]}.')
if np.logical_not(kept_magmoms[np.argmin(energies)]):
    print('WARNING: values of initial and final magnetic moments of the most stable configuration '
          'significantly differ.')

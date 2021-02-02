"""
automag.2_coll.plot_results
===========================

Script which plots results of magnetic relaxations.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

import numpy as np
import matplotlib.pyplot as plt

# increase matplotlib pyplot font size
plt.rcParams.update({'font.size': 18})

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
lines, maginfos = [], []
with open('Fe12O18_relax.txt', 'r') as f:
    for line in f:
        if line[0] != ' ':
            lines.append(line)
        else:
            maginfos.append(line)

# extract the results
configurations, energies, kept_magmoms = [], [], []
for line, maginfo in zip(lines, maginfos):
    values = line.split()
    configurations.append(values[0])
    energies.append(float(values[-1].split('=')[1]))

    if values[2].split('=')[1] == 'NONCONVERGED':
        print(f'WARNING: relax of {values[0]} did not converge.')
    if values[3].split('=')[1] == 'NONCONVERGED':
        print(f'WARNING: energy calculation of {values[0]} did not converge.')

    initial, final = maginfo.split('final_magmoms=')
    initial = initial[initial.index('[') + 1:initial.index(']')]
    final = final[final.index('[') + 1:final.index(']')]
    initial = np.array(initial.split(', '), dtype=float)
    final = np.array(final.split(), dtype=float)
    if max(abs(initial - final)) < 1:
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

# plot results
plt.scatter(repr_configurations[~kept_magmoms], energies[~kept_magmoms], c='r')
plt.scatter(repr_configurations[kept_magmoms], energies[kept_magmoms], c='b')
plt.xlabel('configurations')
plt.ylabel('free energy TOTEN')
plt.savefig('stability.png')
# plt.show()

print(f'The most stable configuration is {configurations[np.argmin(energies)]}.')
if kept_magmoms[np.argmin(energies)] is False:
    print('WARNING: values of initial and final magnetic moments significantly differ.')

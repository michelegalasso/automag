"""
automag.2_eofstate.2_plot_results
===================================

Script which plots results of convergence tests.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

# still to be modified..

import os
import matplotlib

import matplotlib.pyplot as plt
import numpy as np

from ase.io import read


# matplotlib backend and font size
matplotlib.use("TkAgg")
plt.rcParams.update({"font.size": 14})

### START OF INPUT PART ###

# root dir for the convergence test
calc_dir = "/home/michele/EXCHANGE/Sm4H23/3_eofstate"

# variable to plot ("Vol/atom [A^3]" or "Magnetic moment [mu_B]")
variable = "Magnetic moment [mu_B]"

# number of heavy atoms in the unit cell (ignored if plotting magmom)
n_heavy_atoms = 8

### END OF INPUT PART ###

# consistency check
assert variable in ["Vol/atom [A^3]", "Magnetic moment [mu_B]"]

# initialization
pressure_values = []
variable_values = []

for folder in sorted(os.listdir(calc_dir)):
    if folder == "raw_input":
        continue

    # read OUTCAR
    atoms = read(os.path.join(calc_dir, folder, "OUTCAR"))
    pressure = float(folder)

    pressure_values.append(pressure)
    if variable == "Vol/atom [A^3]":
        variable_values.append(atoms.cell.volume / n_heavy_atoms)
    else:
        variable_values.append(atoms.calc.results["magmoms"][0, 0])

plt.plot(pressure_values, variable_values, marker=".")
plt.xlabel("Pressure [GPa]")
plt.ylabel(variable)
# plt.savefig("eofstate.png", bbox_inches="tight")
plt.tight_layout()
plt.show()

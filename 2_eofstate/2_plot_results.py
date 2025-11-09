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


# matplotlib backend and font size
matplotlib.use("TkAgg")
plt.rcParams.update({"font.size": 14})


def read_energy_from_outcar(filename):
    # initialization
    result = None

    # read energy (sigma -> 0)
    with open(filename) as f:
        for line in f:
            if "energy  without entropy" in line:
                result = float(line.split()[-1])

    # check that the energy value has been read
    assert result is not None

    # return energy value
    return result


### START OF INPUT PART ###

# choose the desired mode: 'encut' or 'kgrid'
mode = "kgrid"

# root dir for the convergence test
calc_dir = "/home/michele/EXCHANGE/Sm4H23/1_conv_tests/kgrid"

### END OF INPUT PART ###

# consistency check
assert mode in ["encut", "kgrid"]

if mode == "encut":
    # initialization
    encut_values = []
    energy_values = []

    for i, folder in enumerate(sorted(os.listdir(calc_dir))):
        if i == 0:
            # fetch number of atoms in the unit cell
            with open(os.path.join(calc_dir, folder, "POSCAR"), "r") as f:
                for j, line in enumerate(f):
                    if j == 6:
                        natoms_list = line.split()
                        natoms = sum([int(item) for item in natoms_list])

        if folder == "raw_input":
            continue

        encut = float(folder)
        energy_atom = read_energy_from_outcar(os.path.join(calc_dir, folder, "OUTCAR")) / natoms

        encut_values.append(encut)
        energy_values.append(energy_atom)

    # encut_values = np.array([300, 350, 400, 450, 500, 550, 600, 650, 700, 750])
    # energy_values = np.array([-208.46069675, -208.78697811, -209.04029607, -209.16239793, -209.27367880, -209.33598966,
    #                           -209.37279351, -209.38261711, -209.38414901, -209.38484983]) / 54

    plt.plot(encut_values, energy_values, marker=".")
    plt.xlabel("ENCUT")
    plt.ylabel("Energy [eV/atom]")
    # plt.savefig("ENCUT.png", bbox_inches="tight")
    plt.tight_layout()
    plt.show()

else:
    # initialization
    kpoints_values = []
    sigma_values = []
    energy_values = []

    for i, folder in enumerate(sorted(os.listdir(calc_dir))):
        if i == 0:
            # fetch number of atoms in the unit cell
            with open(os.path.join(calc_dir, folder, "POSCAR"), "r") as f:
                for j, line in enumerate(f):
                    if j == 6:
                        natoms_list = line.split()
                        natoms = sum([int(item) for item in natoms_list])

        if folder == "raw_input":
            continue

        kpoints_string, sigma_string = folder.split("_")

        kpoints_number = 1
        for char in kpoints_string:
            kpoints_number *= int(char)

        energy_atom = read_energy_from_outcar(os.path.join(calc_dir, folder, "OUTCAR")) / natoms

        kpoints_values.append(kpoints_number)
        sigma_values.append(float("0." + sigma_string[1:]))
        energy_values.append(energy_atom)

    # go to numpy arrays
    kpoints = np.array(kpoints_values)
    sigmas = np.array(sigma_values)
    energies = np.array(energy_values)

    # get unique sigma values
    unique_sigmas = np.unique(sigmas)

    # create the plot
    plt.figure()

    for sigma in unique_sigmas:
        # select data for this sigma
        mask = sigmas == sigma
        k = kpoints[mask]
        e = energies[mask]

        # Sort by kpoints for proper line plotting
        sort_idx = np.argsort(k)
        k_sorted = k[sort_idx]
        e_sorted = e[sort_idx]

        plt.plot(np.cbrt(k_sorted), e_sorted, marker='o', label=f'σ = {sigma}')

    # labels and legend
    plt.xlabel('Mean k-points per direction')
    plt.ylabel('Energy [eV/atom]')
    plt.legend(title='Sigma')
    plt.grid(True)
    # plt.savefig("KGRID.png", bbox_inches="tight")
    plt.tight_layout()
    plt.show()

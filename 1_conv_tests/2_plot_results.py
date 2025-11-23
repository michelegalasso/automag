import os
import matplotlib

import matplotlib.pyplot as plt
import numpy as np

from ase.io import read


# matplotlib backend and font size
matplotlib.use("TkAgg")
plt.rcParams.update({"font.size": 14})


### START OF INPUT PART ###

# choose the desired mode: 'encut' or 'kgrid'
mode = "encut"

# root dir for the convergence test
calc_dir = "/home/michele/EXCHANGE/collinear/Tb4H23/1_conv_tests/encut"

### END OF INPUT PART ###

# consistency check
assert mode in ["encut", "kgrid"]

if mode == "encut":
    # initialization
    encut_values = []
    energy_values = []

    for i, folder in enumerate(sorted(os.listdir(calc_dir))):
        if folder != "raw_input":
            encut = float(folder)
            atoms = read(os.path.join(calc_dir, folder, "OUTCAR"))

            energy_atom = atoms.get_total_energy() / atoms.get_number_of_atoms()

            encut_values.append(encut)
            energy_values.append(energy_atom)

    plt.plot(encut_values, energy_values, marker=".")
    plt.xlabel("ENCUT")
    plt.ylabel("Energy [eV/atom]")
    plt.grid(True)
    # plt.savefig("ENCUT.png", bbox_inches="tight")
    plt.tight_layout()
    plt.show()

else:
    # initialization
    kpoints_values = []
    sigma_values = []
    energy_values = []

    for i, folder in enumerate(sorted(os.listdir(calc_dir))):
        if folder != "raw_input":
            kpoints_string, sigma_string = folder.split("_")
            atoms = read(os.path.join(calc_dir, folder, "OUTCAR"))

            kpoints_number = 1
            for char in kpoints_string:
                kpoints_number *= int(char)

            energy_atom = atoms.get_total_energy() / atoms.get_number_of_atoms()

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

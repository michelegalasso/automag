import os
import numpy as np
import matplotlib.pyplot as plt

from ase.io import read
from pymatgen.io.vasp.inputs import Incar


### START OF INPUT PART ###

# root dir for the convergence test
calc_dir = "/home/michele/EXCHANGE/Tb4H23/4_configurations"

# number of heavy atoms in the unit cell (ignored if plotting magmom)
n_heavy_atoms = 8

### END OF INPUT PART ###


def is_comment(line):
    stripped_line = line.strip()
    if stripped_line.startswith("!") or stripped_line.startswith("#"):
        return True

    return False


def rewrite_incar(path_to_incar):
    # initialization
    new_line = ""

    with open(path_to_incar, "r") as original_file:
        with open(path_to_incar + "2", "w") as new_file:
            for original_line in original_file:
                parts = original_line.split()
                if len(parts) > 0 and parts[-1].endswith("\\") and not is_comment(original_line):
                    new_line += original_line.split("\\")[0]
                    write_flag = False
                else:
                    new_line += original_line
                    write_flag = True

                if write_flag:
                    new_file.write(new_line)
                    new_line = ""


def better_sort(state):
    if state == 'nm':
        return 'aa'
    elif state.startswith('fm'):
        return 'aaa' + state[2:]
    else:
        return state


# initialization
mag_states = []
energies = []
colors = []

for folder in sorted(os.listdir(calc_dir), key=better_sort):
    if folder != "raw_input" and folder != "full_relax":
        # get around Pymatgen bug by rewriting the INCAR as INCAR2
        rewrite_incar(os.path.join(calc_dir, folder, 'INCAR'))

        # read INCAR2 and OUTCAR
        incar = Incar.from_file(os.path.join(calc_dir, folder, 'INCAR2'))
        atoms = read(os.path.join(calc_dir, folder, 'OUTCAR'))
        energy = atoms.get_total_energy()

        # evaluate whether magnetic moments significantly changed
        initial_magmom = np.array(incar.get("MAGMOM"))
        final_magmom = atoms.get_magnetic_moments()

        # the collinear case has not been implemented yet
        if initial_magmom.shape[1] != 3 or final_magmom.shape[1] != 3:
            raise NotImplementedError("The collinear case has not been implemented yet.")

        # compute the norms
        initial_magmom_norms = np.linalg.norm(initial_magmom, axis=1)
        final_magmom_norms = np.linalg.norm(final_magmom, axis=1)

        print("\n")
        for j, item in enumerate(final_magmom):
            if j < 8:
                print(f"{item[0]} {item[1]} {item[2]}")

        # blue means no big variation, red means big variation, orange means angle variation
        color = "tab:blue"
        indices = []

        for i, (initial_norm, final_norm) in enumerate(zip(initial_magmom_norms, final_magmom_norms)):
            # if something initialized as zero goes crazy high we signal a big variation
            if np.isclose(initial_norm, 0.0) and final_norm > 0.5:
                color = "tab:red"

            if initial_norm > 0.5:
                # if something initialized high goes below 0.5 we signal a big variation
                if final_norm < 0.5:
                    color = "tab:red"

                # otherwise we write down the index for later analyzing the angle
                else:
                    indices.append(i)

        # do this only if not flagged as big variation
        if color != "tab:red":
            # go to numpy array
            indices = np.array(indices)

            # get the angle between initial and final magmom vectors in degrees
            dot_products = np.sum(initial_magmom[indices] * final_magmom[indices], axis=1)
            norms = initial_magmom_norms[indices] * final_magmom_norms[indices]
            cos_thetas = np.clip(dot_products / norms, -1.0, 1.0)
            angles_deg = np.degrees(np.arccos(cos_thetas))

            if angles_deg.max() > 5.0:
                color = "tab:orange"

        mag_states.append(folder)
        energies.append(energy)
        colors.append(color)

# numpy arrays
mag_states = np.array(mag_states)
energies = np.array(energies) / n_heavy_atoms

# energy differences
energies = energies - energies.min()

# plot
plt.figure(figsize=(16, 9))
plt.bar(mag_states, energies + 1, bottom=-1, color=colors)

plt.xticks(rotation=90)
# plt.ylim(-0.01, 0.03)
plt.ylabel('Energy [eV/atom]')
plt.show()

import os

import numpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from ase.io import read, write
from prettytable import PrettyTable


### START OF INPUT PART ###

# root dir
calc_dir = "/home/michele/EXCHANGE/Tb4H23/3_configurations/settings002/results"

# magnetic atoms
magnetic_atoms = ["Tb"]

### END OF INPUT PART ###

# matplotlib font size
plt.rcParams.update({"font.size": 14})

def better_sort(state):
    if state == "nm":
        return "aa"
    elif state.startswith("fm"):
        return "aaa" + state[2:]
    else:
        return state


def get_initial_magnetic_moments(path, n_atoms):
    # initialization
    prev_site = 0
    values = []

    with open(path, "r") as f:
        for line in f:
            if "magnetic moment" in line:
                fields = line.split()
                site = int(fields[5])

                if prev_site > site:
                    break
                else:
                    values.append(float(fields[7]))
                    prev_site = site

    remaining_sites = n_atoms - prev_site
    values += remaining_sites * [0.0]
    return np.array(values)

# initialization
mag_states = []
visual_states = []
energies = []
colors = []

for folder in sorted(os.listdir(calc_dir), key=better_sort):
    if os.path.isdir(os.path.join(calc_dir, folder)):
        # read QE output
        atoms = read(os.path.join(calc_dir, folder, "pw.scf.out"))
        energy = atoms.get_total_energy()

        # evaluate whether magnetic moments significantly changed
        initial_magmoms = get_initial_magnetic_moments(
            os.path.join(calc_dir, folder, "pw.scf.out"),
            len(atoms),
        )
        final_magmoms = atoms.get_magnetic_moments()

        # select the magnetic atoms
        magnetic_mask = np.zeros_like(initial_magmoms, dtype=bool)
        for element in magnetic_atoms:
            magnetic_mask = np.logical_or(magnetic_mask, atoms.symbols == element)

        # build visual state
        visual_state = ''
        for magmom in initial_magmoms[magnetic_mask]:
            if magmom > 0:
                visual_state += '+ '
            elif magmom < 0:
                visual_state += '- '
            else:
                visual_state += '0 '

        # blue means no big variation, red means big variation
        color = "tab:blue"
        indices = []

        for initial_magmom, final_magmom in zip(initial_magmoms, final_magmoms):
            # catch the big variations
            if np.abs(initial_magmom - final_magmom) > 0.9:
                color = "tab:red"

        mag_states.append(folder)
        visual_states.append(visual_state[:-1])
        energies.append(energy / len(atoms) * 1000)
        colors.append(color)

# numpy arrays
mag_states = np.array(mag_states)
energies = np.array(energies)

# energy differences
energies = energies - energies.min()

# save settings file
atoms = read(os.path.join(calc_dir, "fm1", "pw.scf.out"))
write("settings.vasp", images=atoms, format="vasp")

# save states file
numerical_states = []
energies_for_file = []
for i, (visual_state, energy, color) in enumerate(zip(visual_states, energies, colors)):
    # skip nm (the first state) and those that changed after scf
    if i > 0 and color == "tab:blue":
        items = visual_state.split()

        numerical_state = []
        for item in items:
            if item == "+":
                numerical_state.append(1)
            elif item == "-":
                numerical_state.append(-1)
            else:
                numerical_state.append(0)

        numerical_states.append(numerical_state)
        energies_for_file.append(energy)

numerical_states = np.array(numerical_states)
numpy.savetxt("states.txt", numerical_states, fmt="%3d")
numpy.savetxt("energies.txt", energies_for_file, fmt="%.8f")

# print table
table = PrettyTable()
table.field_names = ["Name", "Magnetic configuration", "Energy [meV/atom]", "Kept magmoms"]
for mag_state, visual_state, energy, color in zip(mag_states, visual_states, energies, colors):
    if color == "tab:blue":
        kept = "Yes"
    else:
        kept = "No"

    table.add_row([mag_state, visual_state, f"{energy:.2f}", kept])

print(table)

# plot
plt.figure(figsize=(16, 9))
bars = plt.bar(mag_states, energies + 50, bottom=-50, color=colors)
plt.bar_label(
    bars,
    labels=[f"{e:.0f}" for e in energies],
    padding=3,
    fontsize=12,
)

plt.xticks(rotation=90)
# plt.ylim(-0.01, 0.03)
plt.ylabel("Energy [meV/atom]")
# plt.grid()

# legend
plt.legend(
    handles=[
        mpatches.Patch(color="tab:blue", label="unchanged"),
        mpatches.Patch(color="tab:red", label="changed"),
    ]
)

plt.show()
# plt.savefig("configurations.png")

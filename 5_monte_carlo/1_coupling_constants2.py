import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from itertools import product
from pymatgen.core.structure import Structure


### START OF INPUT PART ###

# choose the number of nearest neighbors to take into account
n_of_neighbors = 4

# choose the size of the control group
control_group_size = 0.4

# append coupling constant before launching 2_write_vampire_ucf.py
append_coupling_constants = False

# choose the atomic types to be considered magnetic (default transition metals)
magnetic_atoms = ['Tb']

### END OF INPUT PART ###


def nth_neighbor_pairs(structure: Structure, n: int, tol: float = 1e-3):
    """
    Find n-th neighbor atom pairs in a pymatgen Structure, including periodic replicas,
    avoiding double counting under i <-> j exchange.

    Parameters
    ----------
    structure : pymatgen.core.Structure
        Input structure (e.g. read from sublattice.vasp)
    n : int
        Neighbor shell index (1 = first neighbors, 2 = second neighbors, ...)
    tol : float
        Numerical tolerance for distance comparisons (in Å)

    Returns
    -------
    pairs : list of tuple(int, int)
        List of atom index pairs (i, j). The same (i, j) may appear multiple times
        if multiple periodic replicas satisfy the neighbor condition.
    distance : float
        The n-th neighbor distance (Å)
    """

    lattice = structure.lattice
    frac_coords = np.array(structure.frac_coords)
    natoms = len(structure)

    # search replicas in neighboring cells
    shifts = list(product([-1, 0, 1], repeat=3))

    # collect all distances (unique shells)
    all_distances = []

    for i in range(natoms):
        for j in range(natoms):
            if i == j:
                continue
            for shift in shifts:
                delta_frac = frac_coords[j] + np.array(shift) - frac_coords[i]
                d = lattice.get_cartesian_coords(delta_frac)
                dist = np.linalg.norm(d)
                if dist > tol:
                    all_distances.append(dist)

    # identify neighbor shells
    shells = np.unique(np.round(all_distances, 6))
    shells.sort()

    if n > len(shells):
        raise ValueError("Requested neighbor shell does not exist.")

    target_dist = shells[n - 1]

    # now collect pairs at this distance
    pairs = []

    for i in range(natoms):
        for j in range(i + 1, natoms):  # avoid i <-> j double counting
            for shift in shifts:
                delta_frac = frac_coords[j] + np.array(shift) - frac_coords[i]
                d = lattice.get_cartesian_coords(delta_frac)
                dist = np.linalg.norm(d)
                if abs(dist - target_dist) < tol:
                    pairs.append((i, j))

    return pairs, float(np.round(target_dist, 3))


def system(configurations, shells):
    # construct the A matrix
    matrix = []
    for configuration in configurations:
        equation = [1]
        for shell in shells:
            contribution = 0
            for pair in shell:
                contribution += configuration[pair[0]] * configuration[pair[1]]
            equation.append(contribution)

        matrix.append(equation)

    return np.array(matrix)

# initialize variables
structure = None
states = None
energies = None

# read input files from previous step
for item in os.listdir('../3_configurations'):
    rel_path = os.path.join('../3_configurations', item)
    if os.path.isfile(rel_path):
        if item.startswith('settings') and item.endswith('.vasp'):
            structure = Structure.from_file(rel_path)
        if item.startswith('states') and item.endswith('.txt'):
            states = np.loadtxt(rel_path, dtype=int)
        if item.startswith('energies') and item.endswith('.txt'):
            energies = np.loadtxt(rel_path)

if structure is None:
    raise IOError('No settings file found in ../3_configurations folder.')
if states is None:
    raise IOError('No states file found in ../3_configurations folder.')
if energies is None:
    raise IOError('No energies file found in ../3_configurations folder.')

# find out which atoms are magnetic
for element in structure.composition.elements:
    if 'magnetic_atoms' not in globals():
        element.is_magnetic = element.is_transition_metal
    else:
        if element.name in magnetic_atoms:
            element.is_magnetic = True
        else:
            element.is_magnetic = False

# from meV/atom to total energy of the unit cell in eV/atom
energies = energies * structure.num_sites / 1000

# remove non-magnetic atoms from Structure object
non_magnetic_atoms = [element.symbol for element in structure.composition.elements if not element.is_magnetic]
structure.remove_species(non_magnetic_atoms)

# create fit and control group
fit_group_size = 1 - control_group_size
rng = np.random.default_rng(seed=42)
indices = rng.permutation(len(states))

index = int(fit_group_size * len(states))
fitted_idx = indices[:index]
control_idx = indices[index:]

# get neighbors shell
neighbors_shells = []
distances = []
for k in range(n_of_neighbors):
    pairs, dist = nth_neighbor_pairs(structure, k + 1)
    neighbors_shells.append(pairs)
    distances.append(dist)

# get coupling constants on the fitted group
A = system(states[fitted_idx], neighbors_shells)
values = np.linalg.lstsq(A, energies[fitted_idx], rcond=None)

# evaluate coupling constants on the control group
B = system(states[control_idx], neighbors_shells)
predictions = B @ values[0]
PCC = np.corrcoef(predictions, energies[control_idx])

plt.rcParams.update({'font.size': 13})
plt.locator_params(axis='x', nbins=5)
plt.xlabel(f'Heisenberg model energy [meV/atom]')
plt.ylabel(f'DFT energy [meV/atom]')

# print results
if np.linalg.matrix_rank(A) == len(distances) + 1:
    coupling_constants = values[0][1:] * 1.60218e-19
    print(f'distances between neighbors: {distances}')
    print(f'coupling constants: {np.array2string(coupling_constants, precision=8, separator=", ")}')

    if append_coupling_constants:
        with open('input.py', 'a') as f:
            f.write('\n# LINE ADDED BY THE SCRIPT 1_coupling_constants.py')
            f.write(f'\ndistances_between_neighbors = {distances}\n')

            f.write('\n# LINE ADDED BY THE SCRIPT 1_coupling_constants.py')
            f.write(f'\ncoupling_constants = {np.array2string(coupling_constants, precision=8, separator=", ")}\n')

    plt.scatter(
        predictions,
        energies[control_idx],
        label=f'PCC: {PCC[0, 1]:.2f}',
    )
    ax = plt.gca()
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
    plt.legend()
    plt.show()
    # plt.savefig('wrong.png', bbox_inches='tight')
    print(f'PCC: {PCC[0, 1]:.2f}')
else:
    print(f'ERROR: SYSTEM OF {np.linalg.matrix_rank(A)} INDEPENDENT EQUATION(S) '
          f'IN {len(distances) + 1} UNKNOWNS!')

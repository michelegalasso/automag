# name of the poscar file to use in the automag/geometries folder
poscar_file = 'Fe2O3-alpha_primitive.vasp'

# define the atom to introduce as dummy atom
dummy_atom = 'Zn'

# define the position of the dummy atom (from 0 to N_ATOMS-1)
dummy_position = 0

# define the perturbations in eV to apply to the dummy atom
perturbations = [-0.08, -0.05, -0.02, 0.02, 0.05, 0.08]

# define the VASP parameters
params = {
    'xc': 'PBE',
    'setups': 'recommended',
    'prec': 'Accurate',
    'ncore': 4,
    'encut': 820,
    'ediff': 1e-6,
    'ismear': 0,
    'sigma': 0.05,
    'kpts': 20,
    'lmaxmix': 4,
    'nelm': 200,
}

# choose the atomic types to be considered magnetic (default transition metals)
# magnetic_atoms = ['Fe', 'Co']

# choose the magnetic configuration to use for U calculation (default FM-HS)
# configuration = 6 * [4.0] + 6 * [-4.0] + 18 * [0.0]

# name of the poscar file to use in the automag/geometries folder
poscar_file = 'Fe2O3-alpha.vasp'

# define the VASP parameters
params = {
    'xc': 'PBE',
    'setups': 'recommended',
    'prec': 'Accurate',
    'ncore': 4,
    'encut': 820,
    'ediff': 1e-6,  # it was 1e-8
    'ismear': -5,
    'sigma': 0.05,
    'nelm': 200,
    'kpts': 20,
    'lmaxmix': 4,
    'lcharg': False,
    'lwave': False,
    'isym': 0,
    'ldau': True,
    'ldautype': 2,
    'ldaul': [2, -1],
    'ldauu': [4, 0],
    'ldauj': [0, 0],
    'ldauprint': 2,
}

# choose the atomic types to be considered magnetic (default transition metals)
# magnetic_atoms = ['Fe', 'Co']

# choose the magnetic configuration to use for U calculation (default FM-HS)
# configuration = 6 * [4.0] + 6 * [-4.0] + 18 * [0.0]

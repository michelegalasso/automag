# name of the poscar file to use in the automag/geometries folder
poscar_file = 'Fe2O3-alpha_conventional.vasp'

# maximum supercell size for generating distinct magnetic configurations
supercell_size = 1

# choose the absolute value given to up and down spins
high_spin_value = 4

# if included all Wyckoff positions could be either HS or LS
# low_spin_value = 1

# define the VASP parameters
params = {
    'xc': 'PBE',
    'setups': 'recommended',
    'prec': 'Accurate',
    'ncore': 4,
    'encut': 820,
    'ediff': 1e-6,
    'ediffg': -2e-3,
    'ibrion': 2,
    'isif': 3,
    'nsw': 300,
    'potim': 0.2,
    'ismear': 1,
    'sigma': 0.1,
    'nelm': 80,
    'kpts': 20,
    'lmaxmix': 4,
    'lcharg': False,
    'lwave': False,
    'isym': 0,
    'ldau': True,
    'ldautype': 2,
    'ldaul': [2, -1],
    'ldauu': [3.75, 0],
    'ldauj': [0, 0],
    'ldauprint': 2,
}

# choose the atomic types to be considered magnetic (default transition metals)
# magnetic_atoms = ['Mn']

# specify a cutoff for picking only high-spin configurations from output
# lower_cutoff = 1.7

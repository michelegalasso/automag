# name of the poscar file to use in the automag/geometries folder
poscar_file = 'Co_conventional.vasp'

# maximum supercell size for generating distinct magnetic configurations
supercell_size = 2

# choose the absolute value given to up and down spins
spin_value = 2

# define the VASP parameters
params = {
    'xc': 'PBE',
    'setups': 'recommended',
    'prec': 'Accurate',
    'ncore': 4,
    'encut': 470,
    'ediff': 1e-6,
    'ismear': 1,
    'sigma': 0.1,
    'nelm': 200,
    'kpts': 40,
    # 'lmaxmix': 4,
    'lcharg': False,
    'lwave': False,
    'isym': 0,
    # 'ldau': True,
    # 'ldautype': 2,
    # 'ldaul': [2, -1],
    # 'ldauu': [3.75, 0],
    # 'ldauj': [0, 0],
    # 'ldauprint': 2,
}

# choose the atomic types to be considered magnetic (default transition metals)
# magnetic_atoms = ['Mn']

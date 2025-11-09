# name of the poscar file to use in the automag/geometries folder
poscar_file = 'Ni3TeO6_setting005_1x2x1.vasp'

# maximum supercell size for generating distinct magnetic configurations
supercell_size = 1

# choose the absolute values given to up and down spins
spin_values = {
    'Ni': [2],
}

# define the VASP parameters
params = {
    'xc': 'PBE',
    'setups': 'recommended',
    'prec': 'Accurate',
    'ncore': 4,
    'encut': 830,
    'ediff': 1e-6,
    'ismear': 1,
    'sigma': 0.1,
    'nelm': 200,
    'kpts': 20,
    'lmaxmix': 4,
    'lcharg': False,
    'lwave': False,
    'isym': 0,
    'ldau': True,
    'ldautype': 2,
    'ldaul': [2, -1, -1],
    'ldauu': [5.17, 0, 0],
    'ldauj': [0, 0, 0],
    'ldauprint': 2,
}

# specify a cutoff for picking only high-spin configurations from output
# lower_cutoff = 1.7

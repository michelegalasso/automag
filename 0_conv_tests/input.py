# choose the desired mode: 'encut' or 'kgrid'
mode = 'encut'

# name of the poscar file to use in the automag/geometries folder
poscar_file = 'Fe-bcc.vasp'

# define the VASP parameters
params = {
    'xc': 'PBE',
    'setups': 'recommended',
    'istart': 0,
    'prec': 'Normal',
    'ncore': 4,
    'ediff': 1e-6,
    'ismear': 1,
    'lcharg': False,
    'lwave': False,
    'sigma': 0.1,
    'kpts': 30,
    # 'encut': 670,
}

# choose the atomic types to be considered magnetic (default transition metals)
# magnetic_atoms = ['Fe', 'Co']

# choose the magnetic configuration to use for convergence tests (default FM-HS)
# configuration = 6 * [4.0] + 6 * [-4.0] + 18 * [0.0]

# choose the trial values for ENCUT (default from 400 to 1000 eV at steps of 10 eV)
encut_values = range(250, 510, 10)

# choose the trial values for ENCUT (default from 400 to 1000 eV at steps of 10 eV)
# sigma_values = [0.05, 0.1, 0.2]

# choose the trial values for ENCUT (default from 400 to 1000 eV at steps of 10 eV)
# kpts_values = range(20, 110, 10)

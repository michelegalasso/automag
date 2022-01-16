"""
automag.1_lin_response.submit
=============================

Script which submits linear response U calculations.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

from copy import copy
from ase.io import read, write

from common.SubmitFirework import SubmitFirework

# location of the poscar file with the input structure
atoms = read('../input/Fe2O3-alpha.vasp')

# choose the magnetic state to use for U calculation
magmoms = 12 * [4.0] + 18 * [0.0]

# choose the desired mode: 'BARE', 'NSC' or 'SC'
MODE = 'BARE'

# launchdir of the bare run
BARE_DIR = '/cephfs/home/mgalasso/fw_calcs/block_2021-01-09-08-46-32-079078/launcher_2021-01-19-12-25-04-409802'

# change first Fe atom to Ni and save in a poscar file
poscar_file = 'onefakeatom.vasp'
ch_symbols = atoms.get_chemical_symbols()
ch_symbols[0] = 'Ni'
atoms.set_chemical_symbols(ch_symbols)
write(poscar_file, atoms)

# define the VASP parameters
bare_params = {
    'xc': 'PBE',
    'setups': 'recommended',
    'prec': 'Accurate',
    'ncore': 4,
    'encut': 670,
    'ediff': 1e-6,
    'ismear': 0,
    'sigma': 0.2,
    # 'nbands': 148,
    # 'pstress': 1500,
    'kpts': 30,
    'lmaxmix': 4,
}

# submit bare run
if MODE == 'BARE':
    bare_run = SubmitFirework(poscar_file, mode='singlepoint', fix_params=bare_params, magmoms=magmoms)
    bare_run.submit()

# submit non-selfconsistent response
nsc_params = copy(bare_params)
nsc_params['ldau'] = True
nsc_params['ldautype'] = 3
nsc_params['ldaul'] = [2, -1, -1]
nsc_params['icharg'] = 11

perturbations = [-0.08, -0.05, -0.02, 0.02, 0.05, 0.08]
if MODE == 'NSC':
    nsc_run = SubmitFirework(poscar_file, mode='perturbations', fix_params=nsc_params, magmoms=magmoms,
                             pert_values=perturbations, bare_dir=BARE_DIR)
    nsc_run.submit()

# submit selfconsistent response
sc_params = copy(nsc_params)
del sc_params['icharg']

if MODE == 'SC':
    nsc_run = SubmitFirework(poscar_file, mode='perturbations', fix_params=sc_params, magmoms=magmoms,
                             pert_values=perturbations)
    nsc_run.submit()

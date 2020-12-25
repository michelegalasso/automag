"""
automag.1_lin_response.submit
==================================

Script which submits linear response U calculations.

.. codeauthor:: Michele Galasso <michele.galasso@skoltech.ru>
"""

from copy import copy
from ase.io import read, write

from common.SubmitFirework import SubmitFirework

# location of the poscar file with the input structure
atoms = read('../input/Fe2O3-alpha.vasp')

# what are we launching
# MODE = 'BARE'
# MODE = 'NSC'
MODE = 'SC'

# launchdir of the bare run
BARE_DIR = '/cephfs/home/mgalasso/Zfw_calcs/block_2020-12-14-08-02-00-715229/launcher_2020-12-25-08-35-18-089220'

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
    'nbands': 148,
    'kpts': 30,
    'lmaxmix': 4,
}

# submit bare run
if MODE == 'BARE':
    bare_run = SubmitFirework(poscar_file, mode='singlerun', fix_params=bare_params)
    bare_run.submit()

# submit non-selfconsistent response
nsc_params = copy(bare_params)
nsc_params['ldau'] = True
nsc_params['ldautype'] = 3
nsc_params['ldaul'] = [2, -1, -1]
nsc_params['icharg'] = 11

perturbations = [-0.08, -0.05, -0.02, 0.02, 0.05, 0.08]
if MODE == 'NSC':
    nsc_run = SubmitFirework(poscar_file, mode='perturbations', fix_params=nsc_params, pert_values=perturbations,
                             bare_dir=BARE_DIR)
    nsc_run.submit()

# submit selfconsistent response
sc_params = copy(nsc_params)
del sc_params['icharg']

if MODE == 'SC':
    nsc_run = SubmitFirework(poscar_file, mode='perturbations', fix_params=sc_params, pert_values=perturbations)
    nsc_run.submit()
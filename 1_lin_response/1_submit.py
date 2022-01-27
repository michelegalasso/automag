"""
automag.1_lin_response.submit
=============================

Script which submits linear response U calculations.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

from input import *

from ase.io import read, write
from pymatgen.core.structure import Structure

from common.SubmitFirework import SubmitFirework

# full path to poscar file
path_to_poscar = '../geometries/' + poscar_file

# create ase.Atoms object
atoms = read(path_to_poscar)

if 'configuration' not in globals():
    configuration = []
    structure = Structure.from_file(path_to_poscar)
    for atom in structure.species:
        if 'magnetic_atoms' not in globals():
            if atom.is_transition_metal:
                configuration.append(4.0)
            else:
                configuration.append(0.0)
        else:
            if atom in magnetic_atoms:
                configuration.append(4.0)
            else:
                configuration.append(0.0)

# insert dummy atom and save in a poscar file
poscar_file = 'onedummyatom.vasp'
ch_symbols = atoms.get_chemical_symbols()
ch_symbols[dummy_position] = dummy_atom
atoms.set_chemical_symbols(ch_symbols)
write(poscar_file, atoms)

# define the VASP parameters
bare_params = {
    'xc': 'PBE',
    'setups': 'recommended',
    'prec': 'Accurate',
    'ncore': 4,
    'encut': 820,
    'ediff': 1e-6,
    'ismear': 0,
    'sigma': 0.05,
    # 'nbands': 148,
    # 'pstress': 1500,
    'kpts': 20,
    'lmaxmix': 4,
}

# submit calculations
run = SubmitFirework(poscar_file, mode='perturbations', fix_params=bare_params, pert_values=[-0.08, -0.05],
                              magmoms=configuration)
run.submit()

# submit non-selfconsistent response
# nsc_params = copy(bare_params)
# nsc_params['ldau'] = True
# nsc_params['ldautype'] = 3
# nsc_params['ldaul'] = [2, -1, -1]
# nsc_params['icharg'] = 11
#
# perturbations = [-0.08, -0.05, -0.02, 0.02, 0.05, 0.08]
# if MODE == 'NSC':
#     nsc_run = SubmitFirework(poscar_file, mode='perturbations', fix_params=nsc_params, magmoms=magmoms,
#                              pert_values=perturbations, bare_dir=BARE_DIR)
#     nsc_run.submit()
#
# # submit selfconsistent response
# sc_params = copy(nsc_params)
# del sc_params['icharg']
#
# if MODE == 'SC':
#     nsc_run = SubmitFirework(poscar_file, mode='perturbations', fix_params=sc_params, magmoms=magmoms,
#                              pert_values=perturbations)
#     nsc_run.submit()

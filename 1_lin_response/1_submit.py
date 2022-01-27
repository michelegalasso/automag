"""
automag.1_lin_response.submit
=============================

Script which submits linear response U calculations.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

from input import *

from ase.io import read
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
    'kpts': 20,
    'lmaxmix': 4,
}

# submit calculations
run = SubmitFirework(path_to_poscar, mode='perturbations', fix_params=bare_params, pert_values=perturbations,
                     magmoms=configuration, dummy_atom=dummy_atom, dummy_position=dummy_position)
run.submit()

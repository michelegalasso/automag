"""
automag.0_conv_tests.1_submit
=============================

Script which submits convergence tests.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

from input import *

from pymatgen.core.structure import Structure

from common.SubmitFirework import SubmitFirework

# full path to poscar file
path_to_poscar = '../geometries/' + poscar_file

# magnetic configuration to use for the convergence test
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

# convergence test w.r.t. encut
if mode == 'encut':
    if 'encut_values' not in globals():
        encut_values = range(500, 1010, 10)

    convtest = SubmitFirework(path_to_poscar, mode='encut', fix_params=params, magmoms=configuration,
                              encut_values=encut_values)
    convtest.submit()

# convergence test w.r.t. sigma and kpts
if mode == 'kgrid':
    if 'sigma_values' not in globals():
        sigma_values = [item / 100 for item in range(5, 25, 5)]

    if 'kpts_values' not in globals():
        kpts_values = range(20, 110, 10)

    convtest = SubmitFirework(path_to_poscar, mode='kgrid', fix_params=params, magmoms=configuration,
                              sigma_values=sigma_values, kpts_values=kpts_values)
    convtest.submit()

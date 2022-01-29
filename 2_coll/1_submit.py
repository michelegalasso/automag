"""
automag.2_coll.1_submit
=======================

Script which submits collinear single point energy calculations.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

from input import *

import os
import numpy as np

from pymatgen.core.composition import Composition

from common.SubmitFirework import SubmitFirework
from common.insert_elements_in_poscar import insert_elements_in_poscar

COMPOSITION = Composition('Fe12O18')
MAGNETIC_ATOM = 'Fe'
FM_INIT = 12 * [4.0] + 18 * [0.0]
AFM_INIT = 6 * [4.0] + 6 * [-4.0] + 18 * [0.0]

# name of tmp file
TMP_FILENAME = 'tmp.vasp'

filenames = []
for filename in os.listdir('enumlib'):
    if filename.startswith('vasp'):
        filenames.append(filename)
filenames.sort(key=lambda x: int(x.split('.')[1]))

for filename in filenames:
    with open(os.path.join('enumlib', filename), 'r') as f:
        poscar_string = f.read()

    if sum([int(item) for item in poscar_string.split('\n')[5].split()]) != COMPOSITION.num_atoms:
        COMPOSITION *= 2
        FM_INIT = np.repeat(FM_INIT, 2)
        AFM_INIT = np.repeat(AFM_INIT, 2)

    poscar_string = insert_elements_in_poscar(poscar_string, COMPOSITION, MAGNETIC_ATOM)
    with open(TMP_FILENAME, 'w') as f:
        f.write(poscar_string)

    if filename.split('.')[1] == '1':
        # FM configuration
        singlepoint_run = SubmitFirework(TMP_FILENAME, mode='singlepoint', fix_params=params, magmoms=FM_INIT,
                                         configuration='fm')
        singlepoint_run.submit()

    configuration = 'afm' + filename.split('.')[1]
    # AFM configuration
    singlepoint_run = SubmitFirework(TMP_FILENAME, mode='singlepoint', fix_params=params, magmoms=AFM_INIT,
                                     configuration=configuration)
    singlepoint_run.submit()

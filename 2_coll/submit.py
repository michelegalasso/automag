"""
automag.2_coll.submit
=====================

Script which submits collinear relaxations.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

import os

from common.SubmitFirework import SubmitFirework
from pymatgen.core.composition import Composition

COMPOSITION = Composition('Fe12O18')
MAGNETIC_ATOM = 'Fe'

# name of tmp file
TMP_FILENAME = 'tmp.vasp'

params = {
    'xc': 'PBE',
    'setups': 'recommended',
    'prec': 'Accurate',
    'ncore': 4,
    'encut': 670,
    'ediff': 1e-8,
    'ediffg': -1e-3,
    'ibrion': 2,
    'isif': 3,
    'nsw': 300,
    'potim': 0.2,
    'ismear': 1,
    'sigma': 0.1,
    # 'nbands': 148,
    # 'pstress': 1500,
    'kpts': 30,
    'lmaxmix': 4,
    'lcharg': False,
    'lwave': False,
    'isym': 0,
    'ldau': True,
    'ldautype': 2,
    'ldaul': [2, -1],
    'ldauu': [4, 0],
    'ldauj': [0, 0],
    'ldauprint': 2,
}

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

    to_replace = '  '
    replace_with_1 = ''
    replace_with_2 = ''
    for element, abundance in COMPOSITION.items():
        if element.name == MAGNETIC_ATOM:
            to_replace += f'{int(abundance / 2)}   {int(abundance / 2)}   '
        else:
            to_replace += f'{int(abundance)}   '

        replace_with_1 += f'{element.name:>3s} '
        replace_with_2 += f'{int(abundance):3d} '

    poscar_string = poscar_string.replace(to_replace, replace_with_1 + '\n' + replace_with_2)
    with open(TMP_FILENAME, 'w') as f:
        f.write(poscar_string)

    if filename.split('.')[1] == '1':
        magmoms = 30 * [0.0]                        # NM configuration
        relax_run = SubmitFirework(TMP_FILENAME, mode='relax', fix_params=params, magmoms=magmoms,
                                   configuration='nm')
        relax_run.submit()

        magmoms = 12 * [4.0] + 18 * [0.0]           # FM configuration
        relax_run = SubmitFirework(TMP_FILENAME, mode='relax', fix_params=params, magmoms=magmoms,
                                   configuration='fm')
        relax_run.submit()

    configuration = 'afm' + filename.split('.')[1]
    magmoms = 6 * [4.0] + 6 * [-4.0] + 18 * [0.0]   # AFM configuration
    relax_run = SubmitFirework(TMP_FILENAME, mode='relax', fix_params=params, magmoms=magmoms,
                               configuration=configuration)
    relax_run.submit()

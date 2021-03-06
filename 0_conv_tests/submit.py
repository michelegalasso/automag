"""
automag.0_conv_tests.submit
===========================

Script which submits convergence tests.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

from copy import copy

from common.SubmitFirework import SubmitFirework

# choose the desired mode: 'encut' or 'kgrid'
MODE = 'encut'

# location of the poscar file with the input structure
poscar_file = '../input/Fe2O3-alpha.vasp'

# choose the magnetic state to use for convergence tests
magmoms = 12 * [4.0] + 18 * [0.0]

# define the VASP parameters
params = {
    'xc': 'PBE',
    'setups': 'recommended',
    'istart': 0,
    'prec': 'Normal',
    'ncore': 4,
    'ediff': 1e-6,
    'ismear': 1,
    # 'nbands': 148,
    # 'pstress': 1500,
    'lcharg': False,
    'lwave': False,
}

# convergence test w.r.t. encut
if MODE == 'encut':
    params1 = copy(params)
    params1['sigma'] = 0.1
    params1['kpts'] = 30

    encut_values = range(500, 1010, 10)

    convtest = SubmitFirework(poscar_file, mode='encut', fix_params=params1, magmoms=magmoms, encut_values=encut_values)
    convtest.submit()

# convergence test w.r.t. sigma and kpts
if MODE == 'kgrid':
    params2 = copy(params)
    params2['encut'] = 670

    sigma_values = [0.05, 0.1, 0.2]
    kpts_values = range(20, 110, 10)

    convtest = SubmitFirework(poscar_file, mode='kgrid', fix_params=params2, magmoms=magmoms, sigma_values=sigma_values,
                              kpts_values=kpts_values)
    convtest.submit()

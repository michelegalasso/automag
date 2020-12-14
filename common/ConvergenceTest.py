"""
automag.ConvergenceTest
==================================

Class for convergence tests.

.. codeauthor:: Michele Galasso <michele.galasso@skoltech.ru>
"""

import numpy as np

from copy import copy
from ase.io import read
from itertools import product
from fireworks import LaunchPad, Firework, Workflow

from utilities import atoms_to_encode, VaspCalculationTask, WriteOutputTask


launchpad = LaunchPad.from_file('/home/michele/.fireworks/my_launchpad.yaml')


class ConvergenceTest(object):
    def __init__(self, fix_params, encut_values = None, sigma_values = None, kpts_values = None):
        if encut_values is not None:
            assert sigma_values is None
            assert kpts_values is None
            self.var_params = product(encut_values)
        else:
            assert sigma_values is not None
            assert kpts_values is not None
            self.var_params = product(sigma_values, kpts_values)

        self.fix_params = fix_params

    def submit(self):
        params = copy(self.fix_params)
        for values in self.var_params:
            if len(values) == 1:
                mode = 'encut'
                params['encut'] = values[0]
                name = mode + str(values[0])
            elif len(values) == 2:
                mode = 'kgrid'
                params['sigma'] = values[0]
                params['kpts'] = values[1]
                name = mode + str(values[0]) + '-' + str(values[1])
            else:
                raise ValueError('Convergence tests for three parameters simultaneously are not supported.')

            # create an atoms object and encode it
            atoms = read('input/Fe2O3-alpha.vasp')
            magmoms = np.array(12 * [4.0] + 18 * [0.0])         # use the FM configuration
            encode = atoms_to_encode(atoms)

            # single point calculation
            firetask1 = VaspCalculationTask(
                calc_params=params,
                encode=encode,
                magmoms=magmoms,
            )
            sp_firework = Firework(
                [firetask1],
                name=name,
                spec={'_pass_job_info': True, '_queueadapter': {'walltime': '24:00:00'}},
                fw_id=1,
            )

            # write output
            firetask2 = WriteOutputTask(
                system=name,
                filename=f'Zfw_calcs/{atoms.get_chemical_formula()}_{mode}.txt',
                initial_magmoms=magmoms,
                read_enthalpy=False,
            )
            output_firework = Firework(
                [firetask2],
                name='write_output',
                spec={'_queueadapter': {'ntasks': 1, 'walltime': '00:30:00'}},
                fw_id=2
            )

            # package the fireworks into a workflow and submit to the launchpad
            workflow = Workflow([sp_firework, output_firework], name=name, links_dict={1: [2], 2: []})
            launchpad.add_wf(workflow)

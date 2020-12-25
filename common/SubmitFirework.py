"""
automag.common.SubmitFirework
==================================

Class which takes care of submitting FireWorks to the remote database.

.. codeauthor:: Michele Galasso <michele.galasso@skoltech.ru>
"""

import numpy as np

from copy import copy
from ase.io import read
from typing import Union
from itertools import product
from fireworks import LaunchPad, Firework, Workflow

from utilities import atoms_to_encode, VaspCalculationTask, WriteOutputTask


launchpad = LaunchPad.from_file('/home/michele/.fireworks/my_launchpad.yaml')


class SubmitFirework(object):
    def __init__(self, poscar_file: str, mode: str, fix_params: dict, encut_values: Union[list, range] = None,
                 sigma_values: Union[list, range] = None, kpts_values: Union[list, range] = None,
                 pert_values: Union[list, range] = None, bare_dir: str = None):
        if mode == 'encut':
            assert encut_values is not None
            assert sigma_values is None
            assert kpts_values is None
            assert pert_values is None
            self.var_params = product(encut_values)
        elif mode == 'kgrid':
            assert encut_values is None
            assert sigma_values is not None
            assert kpts_values is not None
            assert pert_values is None
            self.var_params = product(sigma_values, kpts_values)
        elif mode == 'perturbations':
            assert encut_values is None
            assert sigma_values is None
            assert kpts_values is None
            assert pert_values is not None
            self.var_params = product(pert_values)
        elif mode == 'singlerun':
            assert encut_values is None
            assert sigma_values is None
            assert kpts_values is None
            assert pert_values is None
            self.var_params = []
        else:
            raise ValueError(f'Value of mode = {mode} not understood.')

        self.mode = mode
        self.bare_dir = bare_dir
        self.fix_params = fix_params
        self.poscar_file = poscar_file

    def submit(self):
        params = copy(self.fix_params)
        if self.var_params:
            for values in self.var_params:
                if len(values) == 1:
                    if self.mode == 'encut':
                        params['encut'] = values[0]
                    else:
                        params['ldauu'] = [values[0], 0, 0]
                        params['ldauj'] = [values[0], 0, 0]
                    name = self.mode + str(values[0])
                elif len(values) == 2:
                    params['sigma'] = values[0]
                    params['kpts'] = values[1]
                    name = self.mode + str(values[0]) + '-' + str(values[1])
                else:
                    raise ValueError('Convergence tests for three parameters simultaneously are not supported.')

                self.add_wflow(params, name, self.mode)
        else:
            self.add_wflow(params, name=self.mode, mode=self.mode)

    def add_wflow(self, params, name, mode):
        # create an atoms object and encode it
        atoms = read(self.poscar_file)
        magmoms = np.array(12 * [4.0] + 18 * [0.0])  # use the FM configuration
        encode = atoms_to_encode(atoms)

        # single point calculation
        if self.bare_dir is not None:
            firetask1 = VaspCalculationTask(
                calc_params=params,
                encode=encode,
                magmoms=magmoms,
                bare_dir=self.bare_dir,
            )
        else:
            firetask1 = VaspCalculationTask(
                calc_params=params,
                encode=encode,
                magmoms=magmoms,
            )

        vasp_firework = Firework(
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
        workflow = Workflow([vasp_firework, output_firework], name=name, links_dict={1: [2], 2: []})
        launchpad.add_wf(workflow)

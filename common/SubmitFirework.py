"""
automag.common.SubmitFirework
=============================

Class which takes care of submitting FireWorks to the remote database.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

import numpy as np

from copy import copy
from ase.io import read
from typing import Union
from itertools import product
from fireworks import LaunchPad, Firework, Workflow

from common.utilities import atoms_to_encode, VaspCalculationTask, WriteOutputTask, WriteChargesTask


# substitute with your launchpad file
launchpad = LaunchPad.from_file('/home/mgalasso/.fireworks/my_launchpad.yaml')


class SubmitFirework(object):
    def __init__(self, poscar_file: str, mode: str, fix_params: dict, magmoms: list,
                 encut_values: Union[list, range] = None, sigma_values: Union[list, range] = None,
                 kpts_values: Union[list, range] = None, pert_values: Union[list, range] = None,
                 configuration: str = None, dummy_atom: str = None, dummy_position: int = None):
        if mode == 'encut':
            assert encut_values is not None
            assert sigma_values is None
            assert kpts_values is None
            assert pert_values is None
            assert dummy_atom is None
            assert dummy_position is None
            self.var_params = product(encut_values)
        elif mode == 'kgrid':
            assert encut_values is None
            assert sigma_values is not None
            assert kpts_values is not None
            assert pert_values is None
            assert dummy_atom is None
            assert dummy_position is None
            self.var_params = product(sigma_values, kpts_values)
        elif mode == 'perturbations':
            assert encut_values is None
            assert sigma_values is None
            assert kpts_values is None
            assert pert_values is not None
            assert dummy_atom is not None
            assert dummy_position is not None
            self.var_params = []
        elif mode in ['relax', 'singlepoint']:
            assert encut_values is None
            assert sigma_values is None
            assert kpts_values is None
            assert pert_values is None
            assert dummy_atom is None
            assert dummy_position is None
            self.var_params = []
        else:
            raise ValueError(f'Value of mode = {mode} not understood.')

        if mode == 'encut':
            self.energy_convergence = True
        else:
            self.energy_convergence = False

        self.magmoms = np.array(magmoms)
        self.mode = mode
        self.fix_params = fix_params
        self.poscar_file = poscar_file
        self.configuration = configuration
        self.pert_values = pert_values
        self.dummy_atom = dummy_atom
        self.dummy_position = dummy_position

    def submit(self):
        params = copy(self.fix_params)
        if self.var_params:
            for values in self.var_params:
                if len(values) == 1:
                    if self.mode == 'encut':
                        params['encut'] = values[0]
                    # else:
                    #     params['ldauu'] = [values[0], 0, 0]
                    #     params['ldauj'] = [values[0], 0, 0]
                    name = self.mode + str(values[0])
                elif len(values) == 2:
                    params['sigma'] = values[0]
                    params['kpts'] = values[1]
                    name = self.mode + str(values[0]) + '-' + str(values[1])
                else:
                    raise ValueError('Convergence tests for three parameters simultaneously are not supported.')

                self.add_wflow(params, name)
        else:
            if self.configuration is not None:
                self.add_wflow(params, self.configuration)
            else:
                self.add_wflow(params, self.mode)

    def add_wflow(self, params, name):
        # create an atoms object and encode it
        if self.mode == 'perturbations':
            atoms = read(self.poscar_file)
            ch_symbols = atoms.get_chemical_symbols()
            atom_ucalc = ch_symbols[self.dummy_position]
            ch_symbols[self.dummy_position] = self.dummy_atom
            atoms.set_chemical_symbols(ch_symbols)
            encode = atoms_to_encode(atoms)
        else:
            atoms = read(self.poscar_file)
            encode = atoms_to_encode(atoms)

        # here we will collect all fireworks of our workflow
        fireworks = []

        if self.mode == 'relax':
            relax_firetask = VaspCalculationTask(
                calc_params=params,
                encode=encode,
                magmoms=self.magmoms,
            )
            relax_firework = Firework(
                [relax_firetask],
                name='relax',
                spec={'_pass_job_info': True},
                fw_id=0
            )
            fireworks.append([relax_firework])

            energy_params = {}
            for key, value in params.items():
                if key not in ['ediffg', 'ibrion', 'isif', 'nsw', 'potim', 'ismear', 'sigma']:
                    energy_params[key] = value
            energy_params['ismear'] = -5
            energy_params['sigma'] = 0.05
            energy_params['nelm'] = 200

            # calculate energy
            if self.magmoms.any():
                sp_firetask = VaspCalculationTask(
                    calc_params=energy_params,
                    magmoms='previous',
                )
            else:
                sp_firetask = VaspCalculationTask(
                    calc_params=energy_params,
                    magmoms=self.magmoms,
                )
        else:
            sp_firetask = VaspCalculationTask(
                calc_params=params,
                encode=encode,
                magmoms=self.magmoms,
            )

        sp_firework = Firework(
            [sp_firetask],
            name='singlepoint',
            spec={'_pass_job_info': True},
            fw_id=1,
        )
        fireworks.append([sp_firework])

        if self.mode == 'perturbations':
            next_id = 2
            nsc_fireworks = []
            sc_fireworks = []
            out_fireworks = []
            for perturbation in self.pert_values:
                nsc_firetask = VaspCalculationTask(
                    calc_params=params,
                    encode=encode,
                    magmoms=self.magmoms,
                    pert_step='NSC',
                    pert_value=perturbation,
                    dummy_atom=self.dummy_atom,
                    atom_ucalc=atom_ucalc,
                )

                nsc_firework = Firework(
                    [nsc_firetask],
                    name='nsc',
                    spec={'_pass_job_info': True},
                    fw_id=next_id,
                )

                nsc_fireworks.append(nsc_firework)
                next_id += 1

                sc_firetask = VaspCalculationTask(
                    calc_params=params,
                    encode=encode,
                    magmoms=self.magmoms,
                    pert_step='SC',
                    pert_value=perturbation,
                    dummy_atom=self.dummy_atom,
                    atom_ucalc=atom_ucalc,
                )

                sc_firework = Firework(
                    [sc_firetask],
                    name='sc',
                    spec={'_pass_job_info': True},
                    fw_id=next_id,
                )

                sc_fireworks.append(sc_firework)
                next_id += 1

                out_firetask = WriteChargesTask(
                    filename='charges.txt',
                    pert_value=perturbation,
                    dummy_atom=self.dummy_atom,
                )
                out_firework = Firework(
                    [out_firetask],
                    name='write_charges',
                    spec={'_queueadapter': {'ntasks': 1, 'walltime': '00:30:00'}},
                    fw_id=next_id,
                )

                out_fireworks.append(out_firework)
                next_id += 1

            fireworks.append(nsc_fireworks)
            fireworks.append(sc_fireworks)
            fireworks.append(out_fireworks)
        else:
            # write output
            output_firetask = WriteOutputTask(
                system=name,
                filename=f"{atoms.get_chemical_formula(mode='reduce')}_{self.mode}.txt",
                initial_magmoms=self.magmoms,
                read_enthalpy=False,
                energy_convergence=self.energy_convergence,
            )
            output_firework = Firework(
                [output_firetask],
                name='write_output',
                spec={'_queueadapter': {'ntasks': 1, 'walltime': '00:30:00'}},
                fw_id=2,
            )
            fireworks.append([output_firework])

        # package the fireworks into a workflow and submit to the launchpad
        flat_fireworks = [fw for sublist in fireworks for fw in sublist]

        links_dict = {}
        for i, level in enumerate(fireworks[:-1]):
            next_level = fireworks[i + 1]
            if len(level) == 1:
                links_dict[level[0].fw_id] = [item.fw_id for item in next_level]
            elif len(next_level) == 1:
                for fw in level:
                    links_dict[fw.fw_id] = [next_level[0].fw_id]
            else:
                for j, fw in enumerate(level):
                    links_dict[fw.fw_id] = [next_level[j].fw_id]

        workflow = Workflow(flat_fireworks, name=name, links_dict=links_dict)
        launchpad.add_wf(workflow)

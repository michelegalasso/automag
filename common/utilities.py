"""
automag.common.utilities
========================

Classes and functions to efficiently handle FireWorks calculations.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

import os
import json
import shutil
import numpy as np

from ase import Atoms
from ase.io import read
from ase.calculators.vasp import Vasp
from fireworks import FiretaskBase, explicit_serialize
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def atoms_to_encode(atoms):
    """
    From Atoms to JSON.

    :param atoms: ase Atoms object.
    :return: JSON string.
    """

    data = {
        'cell': atoms.get_cell().tolist(),
        'scaled_positions': atoms.get_scaled_positions().tolist(),
        'numbers': atoms.get_atomic_numbers().tolist(),
        'pbc': atoms.get_pbc().tolist(),
    }

    return json.dumps(data)


def encode_to_atoms(encode):
    """
    From JSON to Atoms.

    :param encode: JSON string.
    :return: ase Atoms object.
    """
    data = json.loads(encode, encoding='utf-8')

    # construct the Atoms object
    atoms = Atoms(
        cell=data['cell'],
        scaled_positions=data['scaled_positions'],
        numbers=data['numbers'],
        pbc=data['pbc'],
    )

    return atoms


@explicit_serialize
class VaspCalculationTask(FiretaskBase):
    """
    Launches a VASP calculation with the given parameters.
    """
    _fw_name = 'VaspCalculationTask'
    required_params = ['calc_params']
    optional_params = ['encode', 'magmoms', 'bare_dir']

    def run_task(self, fw_spec):
        # if bare_dir is given, copy CHGCAR and WAVECAR from it
        if 'bare_dir' in self:
            shutil.copy(os.path.join(self['bare_dir'], 'CHGCAR'), '.')
            shutil.copy(os.path.join(self['bare_dir'], 'WAVECAR'), '.')

        # if encode is given, use it as input structure
        if 'encode' in self:
            atoms = encode_to_atoms(self['encode'])

        # else, read the input structure from the output of the previous step
        else:
            job_info_array = fw_spec['_job_info']
            prev_job_info = job_info_array[-1]
            atoms = read(os.path.join(prev_job_info['launch_dir'], 'OUTCAR'))

        # if magnetic calculation
        if 'magmoms' in self:
            if self['magmoms'] == 'previous':
                assert 'encode' not in self
                magmoms = atoms.get_magnetic_moments().round()
            else:
                magmoms = np.array(self['magmoms'])

            if magmoms.any():
                # add the necessary VASP parameters
                self['calc_params']['lorbit'] = 11
                self['calc_params']['ispin'] = 2
                atoms.set_initial_magnetic_moments(magmoms)

        # convert any lists from the parameter settings into arrays
        keys = self['calc_params']
        for k, v in keys.items():
            if isinstance(v, list):
                keys[k] = np.asarray(v)

        # initialize and run VASP
        calc = Vasp(**self['calc_params'])
        calc.calculate(atoms)

        # save information about convergence
        with open('is_converged', 'w') as f1:
            for filename in os.listdir('.'):
                if filename.endswith('.out'):
                    with open(filename, 'r') as f2:
                        if calc.converged is False and 'fatal error in bracketing' not in f2.read():
                            f1.write('NONCONVERGED')
                        else:
                            f1.write('converged')


@explicit_serialize
class WriteOutputTask(FiretaskBase):
    """
    Write a simple output file.

    The file contains the USPEX structure ID, information about convergence of each VASP run,
    chemical formula, final energy, initial and final magnetic moments for magnetic calculations.
    """
    _fw_name = 'WriteOutputTask'
    required_params = ['system', 'filename', 'read_enthalpy']
    optional_params = ['initial_magmoms']

    def run_task(self, fw_spec):
        job_info_array = fw_spec['_job_info']
        output_line = '{:12s}  '.format(self['system'])
        output_line += '{:4d}    '.format(job_info_array[0]['fw_id'])

        for job_info in job_info_array:
            with open(os.path.join(job_info['launch_dir'], 'is_converged'), 'r') as f:
                output_line += '{:>6s}={:12s}  '.format(job_info['name'], f.readline())

        atoms_final = read(os.path.join(job_info_array[-1]['launch_dir'], 'OUTCAR'))
        output_line += '{:12s}  '.format(atoms_final.get_chemical_formula(mode='metal'))

        structure = Structure.from_file(os.path.join(job_info_array[-1]['launch_dir'], 'CONTCAR'))
        analyzer = SpacegroupAnalyzer(structure)
        output_line += '{:10s}  '.format(analyzer.get_space_group_symbol())

        if self['read_enthalpy']:
            # get enthalpy from OUTCAR
            with open(os.path.join(job_info_array[-1]['launch_dir'], 'OUTCAR'), 'r') as f:
                for line in f:
                    if 'enthalpy' in line:
                        enthalpy_line = line

            output_line += 'enthalpy={}\n'.format(enthalpy_line.split()[4])
        else:
            output_line += 'energy={}\n'.format(atoms_final.calc.results['free_energy'])

        if 'initial_magmoms' in self:
            magmoms = np.array(self['initial_magmoms'])
            output_line += '          initial_magmoms={}  '.format(np.array2string(magmoms, 1000))

            try:
                magmoms_final = atoms_final.get_magnetic_moments()
            except:
                magmoms_final = np.zeros(len(magmoms))

            output_line += 'final_magmoms={}\n'.format(np.array2string(magmoms_final, 1000))

        with open(os.path.join(os.environ.get('HOME'), self['filename']), 'a') as f:
            f.write(output_line)

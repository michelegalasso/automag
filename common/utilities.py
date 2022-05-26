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
from pymatgen.core.periodic_table import Element
from pymatgen.io.vasp import Incar, Outcar


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
    optional_params = ['encode', 'magmoms', 'pert_step', 'pert_value', 'dummy_atom', 'atom_ucalc']

    def run_task(self, fw_spec):
        # if encode is given, use it as input structure
        if 'encode' in self:
            atoms = encode_to_atoms(self['encode'])

        # else, read the input structure from the output of the previous step
        else:
            job_info_array = fw_spec['_job_info']
            prev_job_info = job_info_array[-1]
            atoms = read(os.path.join(prev_job_info['launch_dir'], 'OUTCAR'))

        if 'pert_step' in self:
            atom_ucalc = Element(self['atom_ucalc'])
            elem_list_sorted, indices = np.unique(atoms.get_chemical_symbols(), return_index=True)
            elem_list = elem_list_sorted[np.argsort(indices)]

            self['calc_params']['ldaul'] = []
            self['calc_params']['ldauu'] = []
            self['calc_params']['ldauj'] = []
            for atom in elem_list:
                if atom == self['dummy_atom']:
                    if atom_ucalc.is_transition_metal:
                        self['calc_params']['ldaul'].append(2)
                    elif atom_ucalc.is_lanthanoid or atom_ucalc.is_actinoid:
                        self['calc_params']['ldaul'].append(3)
                    else:
                        raise ValueError(f"Cannot calculate U for the element {self['atom_ucalc']}")

                    self['calc_params']['ldauu'].append(self['pert_value'])
                    self['calc_params']['ldauj'].append(self['pert_value'])
                else:
                    self['calc_params']['ldaul'].append(-1)
                    self['calc_params']['ldauu'].append(0)
                    self['calc_params']['ldauj'].append(0)

            self['calc_params']['ldau'] = True
            self['calc_params']['ldautype'] = 3

            job_info_array = fw_spec['_job_info']
            prev_job_info = job_info_array[-1]
            shutil.copy(os.path.join(prev_job_info['launch_dir'], 'WAVECAR'), '.')

            # if pert_step is NSC, copy CHGCAR from previous step
            if self['pert_step'] == 'NSC':
                shutil.copy(os.path.join(prev_job_info['launch_dir'], 'CHGCAR'), '.')
                self['calc_params']['icharg'] = 11

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
    required_params = ['system', 'filename', 'read_enthalpy', 'energy_convergence']
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

        if self['energy_convergence']:
            errors = []
            with open(os.path.join(job_info_array[-1]['launch_dir'], 'OUTCAR'), 'r') as f:
                for line in f:
                    if 'kinetic energy error' in line:
                        errors.append(float(line.split()[5]))

            _, indices, counts = np.unique(atoms_final.numbers, return_index=True, return_counts=True)
            num_ions = counts[np.argsort(indices)]
            correction = sum(np.multiply(errors, num_ions))
        else:
            correction = 0

        if self['read_enthalpy']:
            # get enthalpy from OUTCAR
            with open(os.path.join(job_info_array[-1]['launch_dir'], 'OUTCAR'), 'r') as f:
                for line in f:
                    if 'enthalpy' in line:
                        enthalpy_line = line

            output_line += 'enthalpy={}\n'.format(float(enthalpy_line.split()[4]) + correction)
        else:
            output_line += 'energy={}\n'.format(atoms_final.get_potential_energy() + correction)

        if 'initial_magmoms' in self:
            magmoms = np.array(self['initial_magmoms'])
            output_line += '          initial_magmoms={}  '.format(np.array2string(magmoms, 10000))

            try:
                magmoms_final = atoms_final.get_magnetic_moments()
            except:
                magmoms_final = np.zeros(len(magmoms))

            output_line += 'final_magmoms={}\n'.format(np.array2string(magmoms_final, 10000))

        with open(os.path.join(os.environ.get('AUTOMAG_PATH'), 'CalcFold', self['filename']), 'a') as f:
            f.write(output_line)

@explicit_serialize
class WriteChargesTask(FiretaskBase):
    """
    Write charges for U calculation.
    """
    _fw_name = 'WriteChargesTask'
    required_params = ['filename', 'pert_value', 'dummy_atom']

    def run_task(self, fw_spec):
        write_output = True
        job_info_array = fw_spec['_job_info']

        # get dummy index
        with open(os.path.join(job_info_array[-1]['launch_dir'], 'POSCAR'), 'rt') as f:
            poscar_lines = f.readlines()
        elem_list = poscar_lines[0].split()
        num_ions = [int(item) for item in poscar_lines[5].split()]

        dummy_index = 0
        for elem, amount in zip(elem_list, num_ions):
            if elem == self['dummy_atom']:
                if amount == 1:
                    break
                else:
                    raise ValueError('More than one dummy atom in the structure')
            else:
                dummy_index += amount

        # get charges
        charges = []
        incar_bare = Incar.from_file(os.path.join(job_info_array[0]['launch_dir'], 'INCAR'))
        outcar_bare = Outcar(os.path.join(job_info_array[0]['launch_dir'], 'OUTCAR'))
        bare_magmoms = [item['tot'] for item, ref in zip(outcar_bare.magnetization, incar_bare['MAGMOM']) if ref != 0]
        for step in [1, 2]:
            outcar = Outcar(os.path.join(job_info_array[step]['launch_dir'], 'OUTCAR'))
            with open(os.path.join(job_info_array[step]['launch_dir'], 'is_converged'), 'r') as f:
                conv_info = f.readline()

            if 'f' in outcar.charge[dummy_index]:
                charges.append(outcar.charge[dummy_index]['f'])
            else:
                charges.append(outcar.charge[dummy_index]['d'])

            final_magmoms = [item['tot'] for item, ref in zip(outcar.magnetization, incar_bare['MAGMOM']) if ref != 0]
            for bare_magmom, final_magmom in zip(bare_magmoms, final_magmoms):
                if bare_magmom != 0:
                    # if the magnetic moment changes too much, do not write charges in output
                    if final_magmom / bare_magmom < 0.5 or final_magmom / bare_magmom > 2.0 \
                            or conv_info == 'NONCONVERGED':
                        write_output = False

        if write_output:
            with open(os.path.join(os.environ.get('AUTOMAG_PATH'), 'CalcFold', self['filename']), 'a') as f:
                f.write(f"{self['pert_value']:5.2f}  {charges[0]}  {charges[1]}\n")

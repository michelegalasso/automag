"""
automag.common.calculation
==========================

Classes and functions to efficiently handle FireWorks calculations.
"""
import os
import shutil
import getpass
import subprocess


class VASPCalculator(object):
    def __init__(self, specific: str):
        """
        Initializes the class.

        Args:
            specific (str): path to the Specific folder with jobscript and VASP input files
        """
        self.specific = specific
        self.username = getpass.getuser()

    def conv_test(self):
        pass

    def submit(self, parameters: dict, calcfold: str):
        """
        Submits VASP calculation.

        Args:
            parameters (dict): desired INCAR tags (those already present in Specific will be overwritten)
            calcfold (str): path to the calculation folder (will be created if it does not exist)

        Returns:
            str: job id of the submitted calculation
        """
        # create calcfold if it does not exist
        if not os.path.isdir(calcfold):
            os.mkdir(calcfold)

        shutil.copy(os.path.join(self.specific, 'jobscript'), calcfold)
        shutil.copy(os.path.join(self.specific, 'POSCAR'), calcfold)
        shutil.copy(os.path.join(self.specific, 'POTCAR'), calcfold)

        with open(os.path.join(self.specific, 'INCAR'), 'rt') as incar_src:
            with open(os.path.join(calcfold, 'INCAR'), 'wt') as incar_dst:
                for line in incar_src:
                    tag = line.split('=')[0].strip()
                    if tag in parameters:
                        incar_dst.write(' '.join([tag, '=', parameters[tag]]))
                        incar_dst.write('\n')
                    else:
                        incar_dst.write(line)

        if os.path.isfile(os.path.join(self.specific, 'KPOINTS')):
            shutil.copy(os.path.join(self.specific, 'KPOINTS'), calcfold)

        sp = subprocess.run(['sbatch', 'jobscript'], cwd=calcfold, capture_output=True, encoding='utf-8')
        job_id = sp.stdout.split()[-1]

        return job_id

    def has_completed(self, job_id: str):
        """
        Checks if a certain calculation has completed.

        Args:
            job_id: calculation job id

        Returns:
            bool: True if completed, false otherwise
        """
        sp = subprocess.run(['squeue', '-u', self.username], capture_output=True, encoding='utf-8')
        for line in sp.stdout.split('\n'):
            if job_id in line:
                return False

        return True

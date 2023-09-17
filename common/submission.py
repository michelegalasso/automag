"""
automag.common.submission
=========================

Implements local and remote submission through several workload managers.
"""
import os
import shutil
import getpass
import subprocess


class BaseSubmission(object):
    def __init__(self):
        self.username = getpass.getuser()


class Local(BaseSubmission):
    def __init__(self):
        super().__init__()

    @staticmethod
    def submit(command: str, calcfold: str):
        """
        Submits VASP calculation using slurm workload manager.

        Args:
            calcfold (str): path to the calculation folder (will be created if it does not exist)

        Returns:
            str: stdout of the VASP calculation when completed
        """
        sp = subprocess.run(command.split(), cwd=calcfold, capture_output=True, encoding='utf-8')
        return sp.stdout


class Slurm(BaseSubmission):
    def __init__(self):
        super().__init__()

    @staticmethod
    def submit(calcfold: str):
        """
        Submits VASP calculation using slurm workload manager.

        Args:
            calcfold (str): path to the calculation folder (will be created if it does not exist)

        Returns:
            str: job id of the submitted calculation
        """
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

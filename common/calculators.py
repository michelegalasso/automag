"""
automag.common.calculators
==========================

External codes are handled as calculators.
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

    def prepare(self, parameters: dict, calcfold: str):
        """
        Submits VASP calculation.

        Args:
            parameters (dict): desired INCAR tags (those already present in Specific will be overwritten)
            calcfold (str): path to the calculation folder (will be created if it does not exist)

        Returns:
            None
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

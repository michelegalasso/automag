"""
automag.ase.run_vasp
====================

Script which tells ase how to run VASP.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

import os

# vasp 5 local
# load_mkl = 'source /home/mgalasso/intel/compilers_and_libraries_2019.5.281/linux/mkl/bin/mklvars.sh intel64'
# load_mpi = 'source /home/mgalasso/intel/compilers_and_libraries_2019.5.281/linux/mpi/intel64/bin/mpivars.sh intel64'
# exitcode = os.system('{}; {}; mpirun -n 16 /home/mgalasso/softs/vasp.5.4.4/bin/vasp_std'.format(load_mkl, load_mpi))

# vasp 5 module
# exitcode = os.system('module load intel/mkl-11.2.3 mpi/impi-5.0.3 vasp/vasp-5.4.4; mpirun vasp_std')

# vasp 6 module
exitcode = os.system('module load vasp/6.1.1; mpirun vasp_std')

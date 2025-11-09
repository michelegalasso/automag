"""
automag.2_eofstate.1_submit
=============================

Script which submits convergence tests.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

import os
import shutil
import subprocess


### START OF INPUT PART ###

# pressure values in GPa
pressure_values = range(300, 800, 50)

# root dir: it must contain a raw_input folder with INCAR (without ENCUT), POSCAR, POTCAR, KPOINTS and jobscript.sh
calc_dir = "/public/home/ac4fa6jhhb/michele/Tb4H23/3_eofstate"

### END OF INPUT PART ###

files_to_copy = ["INCAR", "POSCAR", "POTCAR", "KPOINTS", "jobscript.sh"]

# consistency check
assert os.path.exists(calc_dir)
assert os.path.exists(os.path.join(calc_dir, "raw_input"))
for filename in files_to_copy:
    assert os.path.exists(os.path.join(calc_dir, "raw_input", filename))

for pressure in pressure_values:
    # clean previous run
    if os.path.exists(os.path.join(calc_dir, str(pressure))):
        shutil.rmtree(os.path.join(calc_dir, str(pressure)))

    os.mkdir(os.path.join(calc_dir, str(pressure)))
    for filename in files_to_copy:
        shutil.copy(os.path.join(calc_dir, "raw_input", filename), os.path.join(calc_dir, str(pressure)))

    # append PSTRESS to INCAR
    with open(os.path.join(calc_dir, str(pressure), "INCAR"), "a") as f:
        f.write("\n\n")
        f.write("# Added by Automag\n")
        f.write("PSTRESS = " + str(pressure * 10) + "\n")

    os.chdir(os.path.join(calc_dir, str(pressure)))
    subprocess.run(["sbatch", "jobscript.sh"])

"""
automag.0_conv_tests.1_submit
=============================

Script which submits convergence tests.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

import os
import shutil
import subprocess


# ENCUT values
encut_values = range(300, 800, 50)

# root dir for the convergence test
calc_dir = "/public/home/ac4fa6jhhb/michele/Tb4H23/1_conv_tests/encut"

# consistency check
assert os.path.exists(calc_dir)
assert os.path.exists(os.path.join(calc_dir, "raw_input"))
for filename in ["INCAR", "POSCAR", "POTCAR", "KPOINTS", "jobscript.sh"]:
    assert os.path.exists(os.path.join(calc_dir, "raw_input", filename))

for encut in encut_values:
    os.mkdir(os.path.join(calc_dir, str(encut)))
    for filename in ["INCAR", "POSCAR", "POTCAR", "KPOINTS", "jobscript.sh"]:
        shutil.copy(os.path.join(calc_dir, "raw_input", filename), os.path.join(calc_dir, str(encut)))

    # append ENCUT to INCAR
    with open(os.path.join(calc_dir, str(encut), "INCAR"), "a") as f:
        f.write("\n\n")
        f.write("# Added by Automag\n")
        f.write("ENCUT = " + str(encut) + "\n")

    subprocess.run(["sbatch jobscript.sh"])

"""
automag.1_lin_response.1_submit
===============================

Script which submits linear response U calculations.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

import os
import shutil
import subprocess


### START OF INPUT PART ###

# root dir (it must contain a raw_input folder with INCAR (without ENCUT), POSCAR, POTCAR, KPOINTS and jobscript.sh)
calc_dir = "/public/home/ac4fa6jhhb/michele/Tb4H23/2_hubbard_u"

# choose mode ("groundstate", "scf" or "nscf")
mode = "scf"

# define the position of the dummy atom (from 0 to N_ATOMS-1, ignored in "groundstate" mode)
dummy_position = 0

# define the perturbations in eV to apply to the dummy atom (ignored in "groundstate" mode)
perturbations = [-0.2, -0.15, -0.1, 0.1, 0.15, 0.2]

### END OF INPUT PART ###

files_to_copy = ["INCAR", "POSCAR", "POTCAR", "KPOINTS", "jobscript.sh"]

# consistency check
assert mode in ["groundstate", "scf", "nscf"]
assert os.path.exists(calc_dir)
assert os.path.exists(os.path.join(calc_dir, "raw_input"))
for filename in files_to_copy:
    assert os.path.exists(os.path.join(calc_dir, "raw_input", filename))

if mode == "groundstate":
    # clean previous run
    if os.path.exists(os.path.join(calc_dir, "groundstate")):
        shutil.rmtree(os.path.join(calc_dir, "groundstate"))

    os.mkdir(os.path.join(calc_dir, "groundstate"))
    for filename in files_to_copy:
        shutil.copy(os.path.join(calc_dir, "raw_input", filename), os.path.join(calc_dir, "groundstate"))

    os.chdir(os.path.join(calc_dir, "groundstate"))
    subprocess.run(["sbatch", "jobscript.sh"])

elif mode == "nscf":
    # clean previous run
    if os.path.exists(os.path.join(calc_dir, "nscf")):
        shutil.rmtree(os.path.join(calc_dir, "nscf"))

    os.mkdir(os.path.join(calc_dir, "nscf"))
    for perturbation in perturbations:
        dir_name = str(perturbation).replace(".", "_").replace("-", "m")
        os.mkdir(os.path.join(calc_dir, "nscf", dir_name))

        for filename in files_to_copy:
            shutil.copy(os.path.join(calc_dir, "raw_input", filename), os.path.join(calc_dir, "nscf", dir_name))

        with open(os.path.join(calc_dir, "nscf", dir_name, "INCAR"), "a") as f:
            f.write("\n\n# Added by Automag\n")
            f.write("LDAU = .TRUE.\n")
            f.write("LDAUTYPE = 3\n")
            f.write("LDAUL = 3 -1 -1\n")
            f.write(f"LDAUU = {perturbation} 0.0 0.0\n")
            f.write(f"LDAUJ = {perturbation} 0.0 0.0\n")
            f.write("ICHARG = 11\n")

        # copy CHGCAR and WAVECAR from the ground state
        shutil.copy(os.path.join(calc_dir, "groundstate", "CHGCAR"), os.path.join(calc_dir, "nscf", dir_name))
        shutil.copy(os.path.join(calc_dir, "groundstate", "WAVECAR"), os.path.join(calc_dir, "nscf", dir_name))

        os.chdir(os.path.join(calc_dir, "nscf", dir_name))
        subprocess.run(["sbatch", "jobscript.sh"])

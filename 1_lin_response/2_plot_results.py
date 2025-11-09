"""
automag.1_lin_response.2_plot_results
=====================================

Script which plots results of linear response U calculation.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

import os
import matplotlib

import matplotlib.pyplot as plt
import numpy as np

from numpy.polynomial import Polynomial


# matplotlib backend and font size
matplotlib.use("TkAgg")
plt.rcParams.update({"font.size": 14})

def read_charge_from_outcar(filename):
    # initialization
    result = None

    # read charge on 1st atom
    flag = False
    with (open(filename) as f):
        for line in f:
            if "aborting loop because EDIFF is reached" in line:
                flag = True
            if flag and len(line.split()) != 0 and line.split()[0] == "1":
                result = float(line.split()[-2])
                break

    # check that the energy value has been read
    assert result is not None

    # return energy value
    return result


### START OF INPUT PART ###

# root dir for the convergence test
calc_dir = "/home/michele/EXCHANGE/Tb4H23/2_hubbard_u"

### END OF INPUT PART ###

# initialization
perturbations = []

for dir_name in os.listdir(os.path.join(calc_dir, "nscf")):
    # get perturbation value
    perturbation = float(dir_name.replace("m", "-").replace("p", "."))
    perturbations.append(perturbation)

# initialization
nscf_responses = np.empty_like(perturbations)
scf_responses = np.empty_like(perturbations)

for step in ["nscf", "scf"]:
    for dir_name in os.listdir(os.path.join(calc_dir, step)):
        perturbation = float(dir_name.replace("m", "-").replace("p", "."))
        response = read_charge_from_outcar(os.path.join(calc_dir, step, dir_name, "OUTCAR"))
        index = perturbations.index(perturbation)

        if step == "nscf":
            nscf_responses[index] = response
        else:
            scf_responses[index] = response

response = read_charge_from_outcar(os.path.join(calc_dir, "groundstate", "OUTCAR"))
perturbations = np.array(perturbations)

# append response for zero perturbation
perturbations = np.append(perturbations, [0.0])
nscf_responses = np.append(nscf_responses, [response])
scf_responses = np.append(scf_responses, [response])

# create figure
plt.figure()

# get indices
indices = np.argsort(perturbations)

# compute angular coefficients
p_nscf = Polynomial.fit(perturbations[indices], nscf_responses[indices], 1)
p_scf = Polynomial.fit(perturbations[indices], scf_responses[indices], 1)
q_nscf, m_nscf = p_nscf.convert()
q_scf, m_scf = p_scf.convert()

# plot
plt.plot(perturbations[indices], nscf_responses[indices], marker = ".", label="nscf")
plt.plot(perturbations[indices], perturbations[indices] * m_nscf + q_nscf, label="nscf interp")
plt.plot(perturbations[indices], scf_responses[indices], marker = ".", label="scf")
plt.plot(perturbations[indices], perturbations[indices] * m_scf + q_scf, label="scf interp")

# labels and legend
plt.xlabel('V [eV]')
plt.ylabel('Number of f-electrons')
plt.legend()
plt.grid(True)
# plt.savefig("UCALC.png", bbox_inches="tight")
plt.tight_layout()
plt.show()

# print result
print(f"U = 1/{m_scf:.4f} - 1/{m_nscf:.4f} = {1/m_scf - 1/m_nscf:.1f}")

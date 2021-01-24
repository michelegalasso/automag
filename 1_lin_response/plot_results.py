"""
automag.1_lin_response.plot_results
===================================

Script which plots results of linear response U calculation.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

import numpy as np
import matplotlib.pyplot as plt

from scipy import stats

# increase matplotlib pyplot font size
plt.rcParams.update({'font.size': 18})

# perturbations used for the VASP runs
perturbations = np.array([-0.08, -0.05, -0.02, 0.02, 0.05])

# number of d electrons on the first Fe site
nscf = np.array([5.857, 5.914, 5.974, 6.057, 6.121])
scf = np.array([5.994, 6.003, 6.010, 6.019, 6.025])

plt.figure(figsize=(16, 9))

slope_nscf, intercept_nscf, r_value_nscf, p_value_nscf, std_err_nscf = stats.linregress(perturbations, nscf)
slope_scf, intercept_scf, r_value_scf, p_value_scf, std_err_scf = stats.linregress(perturbations, scf)

plt.plot(perturbations, nscf, 'ro', label=f'NSCF (slope {slope_nscf:.4f})')
plt.plot(perturbations, scf, 'bo', label=f'SCF (slope {slope_scf:.4f})')

nscf_line = [value * slope_nscf + intercept_nscf for value in perturbations]
scf_line = [value * slope_scf + intercept_scf for value in perturbations]
plt.plot(perturbations, nscf_line, 'r-')
plt.plot(perturbations, scf_line, 'b-')

plt.xlabel('α (eV)')
plt.ylabel('d-electrons on first Fe site')
plt.legend()
# plt.show()
plt.savefig('Ucalc.png')

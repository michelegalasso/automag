"""
automag.1_lin_response.2_plot_results
=====================================

Script which plots results of linear response U calculation.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from scipy import stats

# increase matplotlib pyplot font size
plt.rcParams.update({'font.size': 18})

calcfold = os.path.join(os.environ.get('AUTOMAG_PATH'), 'CalcFold')
data = np.loadtxt(os.path.join(calcfold, 'charges.txt'))

perturbations = data[:, 0]
nscf = data[:, 1]
scf = data[:, 2]

plt.figure(figsize=(16, 9))

slope_nscf, intercept_nscf, r_value_nscf, p_value_nscf, std_err_nscf = stats.linregress(perturbations, nscf)
slope_scf, intercept_scf, r_value_scf, p_value_scf, std_err_scf = stats.linregress(perturbations, scf)

print(f'U = {(1/slope_scf) - (1/slope_nscf):4.2f}')

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

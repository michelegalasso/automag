# choose the configuration to use for Monte Carlo simulation
configuration = 'afm1'

# choose the cutoff radius in Angstrom for neighbor search
cutoff_radius = 3.8

# choose the size of the control group
control_group_size = 0.4

# append coupling constant before launching 2_write_vampire_ucf.py
append_coupling_constants = False

# choose the atomic types to be considered magnetic (default transition metals)
# magnetic_atoms = ['Mn']

# LINE ADDED BY THE SCRIPT 1_coupling_constants.py
distances_between_neighbors = [2.9, 2.97, 3.36, 3.7]

# LINE ADDED BY THE SCRIPT 1_coupling_constants.py
coupling_constants = [-9.78540181e-22, -9.31084685e-22, -7.13999252e-21, -4.57653205e-21]

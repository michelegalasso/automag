# Automag
An automatic workflow software for calculating the ground collinear magnetic
state of a given structure and for estimating the critical temperature of the
magnetically ordered to paramagnetic phase transition.

## Installation
Automag is meant to be run on a computing cluster. The first step in the
installation is to clone the repository

`git clone https://github.com/michelegalasso/automag.git`

I assume that you have Python 3 installed on your system, accessible through
the `python` command, with the `wheel` and the `venv` packages installed.
After cloning, go to the `automag` directory which has appeared and initialize
a Python virtual environment

`cd automag`  
`python -m venv .venv`

After that, activate your virtual environment. If you are working on an Unix-like
system, you can do it with the command

`source .venv/bin/activate`

Finally, you need to install all the necessary Python dependencies for Automag
with the command

`pip install -r requirements.txt`

Before using Automag, make sure that enumlib is installed on your system and
that that your command line has access to the commands `enum.x` and `makeStr.py`.
In addition, Automag needs to know how to use VASP, so you need to edit the file
`automag/ase/run_vasp.py` for Automag to correctly load the MKL and MPI libraries
(if needed) and call the `vasp_std` executable on your system. Then add the
following lines to your `~/.bashrc` file

`export PYTHONPATH=/PATH/TO/automag:$PYTHONPATH`  
`export VASP_SCRIPT=/PATH/TO/automag/ase/run_vasp.py`  
`export VASP_PP_PATH=/PATH/TO/pp`  
`export AUTOMAG_PATH=/PATH/TO/automag`

obviously replacing `/PATH/TO` with the actual path to these files or directories.
The `pp` folder is the VASP pseudopotential library and it should contain two
subfolders named, respectively, `potpaw_LDA` and `potpaw_PBE`.

Before starting to use Automag, you should also make sure to have the FireWorks
library correctly pointing to a MongoDB database and configured for launching
jobs through a queue management system (for more detailed information refer
to the FireWorks documentation). Last but not least, open the file
`automag/common/SubmitFirework.py` and edit line 23 with the location of your
`my_launchpad.yaml` file. Now you are ready to use Automag.

## Convergence tests

When studying a magnetic structure, you may want to start from convergence tests.
Automag allows you to check for convergence of the VASP parameter `ENCUT`, which
is the energy cut-off of the plane wave basis set, and to simultaneously check for
convergence of the two parameters `SIGMA` and `kpts`, which are the electronic
smearing parameter and the k-mesh resolution parameter for Brillouin zone sampling,
respectively.

In order to run convergence tests, go to the folder `0_conv_tests` and insert the
following input parameters in the file `input.py`:

- `mode` can be "encut" for convergence tests with respect to `ENCUT` or "kgrid"
for convergence tests with respect to `SIGMA` and `kpts`;
- `poscar_file` is the name of the file in POSCAR format which contains the input
geometry and which has been put in the folder `automag/geometries`;
- `params` is a collection of VASP parameters to be used during single-point
energy calculations.

In addition, the following optional parameters can be specified:

- `magnetic_atoms` contains the atomic types to be considered magnetic (defaults
to transition metals);
- `configuration` is the magnetic configuration to use for convergence tests
(defaults to ferromagnetic high-spin);
- `encut_values` contains the trial values for `ENCUT` (defaults to the interval
[500, 1000] eV at steps of 10 eV);
- `sigma_values` contains the trial values for `SIGMA` (defaults to the interval
[0.05, 0.20] eV at steps of 0.05 eV);
- `kpts_values` contains the trial values for `kpts` (defaults to the interval
[20, 100] A^-1 at steps of 10 A^-1).

Once you have inserted the input parameters in the file `input.py`, launch the
script `1_submit.py` using the `python` executable of your virtual environment and
a number of workflows containing single-point VASP calculations will be saved in
your remote database. These calculations need to be run in a separate directory
called `CalcFold`. You can enter it and then submit the jobs with the following
commands, being in the `automag` directory

`cd CalcFold`  
`nohup qlaunch -r rapidfire -m 10 --nlaunches=infinite &`

This will enter the directory `CalcFold` and it will invoke the `qlaunch`
command in the background, which constantly checks for new calculations in the
remote database and submits them to the queue management system of your cluster.
The command line option `-m 10` means that `qlaunch` will allow a maximum number
of 10 jobs to be in the queue of your system at the same time. If 10 or more
jobs are waiting or running in your queue, `qlaunch` will not submit more jobs
until their total number becomes less than 10. If you wish, you can change this
parameter to a different number. In the following, I assume that the `qlaunch`
process is  always working in the background.

After all calculations have been completed, you can launch the script named
`2_plot_results.py` in the `0_conv_tests` folder. It will read the output file
that Automag wrote in `CalcFold` and it will produce a plot of the parameters
under study versus energy. In addition, the script will also print on screen the
values of the parameters under study for which the error in energy is less than
1 meV/atom with respect to the most accurate result.

## Calculation of the electronic correlation parameter U by linear response

The linear response formalism for the calculation of the Hubbard U is based on
the application of a series of small perturbations to the first magnetic atom.
During the whole process, the perturbed atom must be treated independently from
the other atoms of the same species, and we achieve this using a simple trick: we
change the chemical identity of the atom to which we want to apply perturbations
to a dummy atomic species, but we place the POTCAR file of the original atomic
species in the `pp` folder corresponding to the dummy atom. For example, if we are
calculating the value of the Hubbard U for Fe in the system Fe12O18 and we choose
Zn as dummy atom, we need to ensure that the same POTCAR file of Fe is used also
for Zn in order to have, instead of Zn, a Fe atom that is treated independently
from the others. This can be achieved with the following commands

`cd $VASP_PP_PATH/potpaw_PBE/Zn`  
`mv POTCAR _POTCAR`  
`cp ../Fe/POTCAR .`

In this way the Automag code, and in particular the ASE library, will treat the
system as if it was ZnFe11O18, but when they will look for the POTCAR file in the
`Zn` folder, they will find the POTCAR of Fe, so the system will remain Fe12O18.
Before launching this calculation, go to the directory `1_lin_response` and set
the necessary parameters in the file `input.py`:

- `poscar_file` is the name of the file in POSCAR format which contains the input
geometry and which has been put in the folder `automag/geometries`;
- `dummy_atom` is the name of the dummy atomic species to use for the atom which
is subject to perturbations (do not forget to manually put the right POTCAR file
in the `pp` folder for this atom);
- `dummy_position` is the position of the dummy atom in your POSCAR file (from 0
to N - 1, where N is the number of atoms in the unit cell);
- `perturbations` contains the values of the perturbations to apply to the chosen
atom in eV;
- `params` is a collection of VASP parameters to be used during single-point
energy calculations.

In addition, the following optional parameters can be specified:

- `magnetic_atoms` contains the atomic types to be considered magnetic (defaults
to transition metals);
- `configuration` is the magnetic configuration to use for U calculation (defaults
to ferromagnetic high-spin).

Once the input parameters have been inserted in the file `input.py`, you can
launch the script `1_submit.py` in order to save the necessary VASP jobs to the
remote database. You will see that the instance of `qlaunch` which is running in
the background will immediately send these jobs to the queue management system
of your cluster. When all calculations are completed, you will find in `CalcFold`
the file `charges.txt`, containing the amount of electrons on the partially
occupied shell of the chosen atom for each value of the applied perturbation,
for both the selfconsistent and the non-selfconsistent runs. Now you can execute
the script `2_plot_results.py`, which will plot the selfconsistent and the
non-selfconsistent responses, it will interpolate them as straight lines to the
least squares and it will calculate their slopes. The value of U is obtained from
U = 1/X - 1/X0, where X and X0 are the selfconsistent and the non-selfconsistent
slopes, respectively.

## Search for the most stable magnetic state

The search for the ground collinear magnetic state consists in generating a
number of trial configurations and in computing their single-point energy, in
order to determine which is the most thermodynamically stable. The trial
configurations differ from each other only by the choice of the unit cell and by
the initialization of the magnetic moments. Note that the value of the magnetic
moment on each atom can change during the single-point energy calculation.
Automag generates trial configurations by separately initializing each Wyckoff
position occupied by magnetic atoms in a ferromagnetic (FM), antiferromagnetic
(AFM) or non magnetic (NM) fashion, taking into account all possible combinations.
A completely non-magnetic (NM) configuration is also generated. In this way,
overall ferrimagnetic (FiM) states are allowed if the magnetic atoms occupy more
than one Wyckoff position. For each magnetic atom, one or two absolute values for
the magnetization can be given in input. In the first case, the given value is
used for initializing all the spin-up and spin-down states in the configurations
generated by Automag. Conversely, if two separate values are given for high-spin
(HS) and low-spin (LS) states, then each Wyckoff position occupied by that
magnetic atom is separately initialized in a HS or LS fashion, taking into account
all possible combinations. In order to run such a calculation, go to the directory
`2_coll` and set the necessary input parameters in the file `input.py`:

- `poscar_file` is the name of the file in POSCAR format which contains the input
geometry and which has been put in the folder `automag/geometries`;
- `supercell_size` is the maximum supercell size for generating distinct magnetic
configurations, in multiples of the input structure;
- `spin_values` a maximum of two values (HS and LS) of the magnetization in Bohr
magnetons for each magnetic atom, used to initialize spin-up and spin-down states;
- `params` is a collection of VASP parameters to be used during single-point
energy calculations.

In addition, the following optional parameter can be specified:

- `lower_cutoff` is the minimum value in Bohr magnetons to which a magnetic moment
can fall in order for the corresponding configuration to be used for estimating
the critical temperature of the material (defaults to zero).

Once the input parameters have been inserted in the file `input.py`, you can
launch the script `1_submit.py` in order to save the necessary VASP jobs to the
remote database. The running instance of `qlaunch` will send these jobs to the
queue management system of your cluster. Once all calculations have been
completed, you can launch the script `2_plot_results.py` which will produce a
number of files containing the histogram plot of the obtained energies for all
trial configurations that successfully completed the single-point energy
calculation. The script also prints on screen the name of the configuration
with lowest energy.

## Calculation of the critical temperature

Automag can calculate the critical temperature of the magnetically ordered to
paramagnetic phase transition from a Monte Carlo simulation of the effective
Hamiltonian which describes the magnetic interaction. Automag computes the
coupling constants of the corresponding Heisenberg model and provides all the
necessary input files to run the Monte Carlo simulation with the VAMPIRE software
package. The simulation itself needs to be run by the user, while Automag can be
used to plot and to fit the results. The accuracy of the Heisenberg model is
evaluated by computing the Pearson Correlation Coefficient (PCC) between the
DFT energies and the predicted energies of a control group of magnetic
configurations. It is worth noting that this approach can be applied only if all
magnetic atoms have the same absolute value of the magnetization. In order to
estimate the critical temperature with Automag, go to the folder `3_monte_carlo`
and set the necessary parameters in the file `input.py`:

- `configuration` is the name of the magnetic configuration to use for the Monte
Carlo simulation (usually you want to put here the name of the configuration with
lowest energy, obtained at the previous step);
- `cutoff_radius` is the maximum distance between two magnetic atoms to be
considered as interacting neighbors;
- `control_group_size` is the relative size of the control group used to evaluate
the accuracy of the Heisenberg model;
- `append_coupling_constants` is a boolean flag which tells Automag whether or not
to append the computed values of the coupling constants to the file `input.py`,
which are needed by the script `2_write_vampire_ucf.py`.

In addition, the following optional parameter can be specified:

- `magnetic_atoms` contains the atomic types to be considered magnetic (defaults
to transition metals).

Once the input parameters have been inserted in the file `input.py`, you can
launch the script `1_coupling_constants.py`. It will compute the coupling
constants for the given cutoff radius and it will print their values on screen.
In addition, it will create a file `model.png`, which contains a plot of the
Heisenberg model energies versus the DFT energies for all the configurations in
the control group. The values of the distances between neighbors, the amounts of
neighboring pairs of magnetic atoms in the unit cell at each distance and the
value of the PCC are also printed on screen.

We suggest to run the script `1_coupling_constants.py` a couple of times with the
`append_coupling_constants` flag set to `False` and with different values of the
`cutoff_radius`, in order to investigate how many neighbors you need to include
for obtaining a well-converged Heisenberg model. Once you are satisfied with the
model's accuracy, run the script a last time with the `append_coupling_constants`
flag set to `True` and Automag will append the values of the coupling constants
and the distances between neighbors to the file `input.py`. Now you are ready to
run the script `2_write_vampire_ucf.py`, which will read the file `input.py` and
will produce a VAMPIRE unit cell file `vamp.ucf`. Now you can run VAMPIRE on your
cluster using the unit cell file written by Automag, an input file and a material
file specific for your problem. Automag contains a sample input file and a sample
material file in the folder `3_monte_carlo/vampire_input`, which can be simply
edited and adapted to the problem under study. Once the VAMPIRE run is done, you
can copy the `output` file in the folder `3_monte_carlo` and run the last script
`3_plot_results.py`. It will produce a plot of the mean magnetization length (of
the spin-up channel for antiferromagnetic materials) versus temperature obtained
from the Monte Carlo simulation and it will fit the data using the analytical
expression of the mean magnetization length, obtaining the values of the critical
temperature and of the critical exponent. The fitted values of these two
parameters are printed on screen.

# Automag
Automatic search for the most stable magnetic state of a
given structure. From the geometry of a structure, you can
find its most stable collinear magnetic configuration and
you can get an estimate of the critical temperature of the
magnetically ordered to paramagnetic phase transition.
At the moment, the code can treat structures with only
one magnetic atomic type in the unit cell.

## Installation
Start by installing Automag on your local machine.
I suppose that you have Python 3 installed in your system,
and that it is accessible through the command `python`.
Give the following commands

1.  `git clone https://github.com/michelegalasso/automag.git`
2.  `cd automag`
3.  `python -m venv .venv`
4.  `.venv\Scripts\activate` or for Linux just use `source .venv/bin/activate`
5.  `pip install -r requirements.txt`

Make sure that enumlib (https://github.com/msg-byu/enumlib)
is installed on your local machine and that your command
line has access to the commands `enum.x` and `makeStr.py`.
Then repeat the above operations on your supercomputer,
in order to install Automag there as well. On your
supercomputer, you need also to create a folder named
`fw_calcs` in your home directory. Here Automag will save
the results of calculations. Automag needs to know how to
launch VASP, so edit the file `automag/ase/run_vasp.py`
and add the following lines to your `~\.bashrc` file

- `export PYTHONPATH=/PATH/TO/automag:$PYTHONPATH`
- `export VASP_SCRIPT=/PATH/TO/automag/ase/run_vasp.py`
- `export VASP_PP_PATH=/PATH/TO/pp`

obviously replacing `/PATH/TO` with the actual path to these
files or directories. The `pp` folder is the VASP
pseudopotential library and it should contain two subfolders
named, respectively, `potpaw_LDA` and `potpaw_PBE`.
Depending on the editor that you use on your local machine
in order to modify and run the code, you may need to add
the first of these three lines to your `~\.bashrc` file
in your local machine as well. With the editor Pycharm
(https://www.jetbrains.com/pycharm) this is not necessary.

Before starting to use Automag, you should also make sure
to have the FireWorks libraries on your supercomputer and
your local machine pointing to the same MongoDB database,
and your FireWorks library on your supercomputer must be
configured for launching jobs through a queue management
system (for more detailed information refer to the
FireWorks documentation). Last but not least, open the file
`automag/common/SubmitFirework.py` and edit line 22 with
the location of your `my_launchpad.yaml` file.

Now you are ready to use Automag.

## Convergence tests

When studying a magnetic structure, you may want to start
from convergence tests. Automag allows you to check for
convergence of the VASP parameter `ENCUT` (the energy cutoff
of the plane wave basis set) and to simultaneously check for
convergence of the two parameters `SIGMA` and `kpts`, which
are the smearing parameter of the electronic temperature
and the k-mesh for Brillouin zone sampling, respectively.

On your local machine, open the file `submit.py` in the
`0_conv_tests` folder. At line 15 choose the desired mode:
'encut' or 'kgrid', at line 18 put the location of the file
which contains your input structure in POSCAR format, at
line 21 you can choose the magnetic state to use for spin
initialization, at lines 25-35 you can tune the other VASP
parameters used during the convergence test and at lines
41-44 and 52-55 you can choose which trial parameters to
use for each convergence test, respectively.

When you are satisfied with your `submit.py` file on your
local machine, launch it in order to save the convergence
test calculations on your remote database. Then go to your
supercomputer and, in the `fw_calcs` folder, give the command

`nohup qlaunch -r rapidfire -m 10 --nlaunches=infinite &`

This will invoke the `qlaunch` command in the background,
which constantly checks for new calculations in the remote
database and submits them to the queue management system.
In the following, I assume that the `qlaunch` process is
always working in the background on your supercomputer.

When calculations are done, copy the output file which you
find in the `fw_calcs` folder in the `0_conv_tests` folder
on your local machine and open the file `plot_results.py`.
At line 17 choose the desired mode and at line 54put the
correct chemical formula of the system under study. Then
launch the script. A PNG file will be created with a visual
representation of the results of the convergence test.

## Calculation of the electronic correlation parameter U by linear response

On your local machine, open the file `submit.py` in the
`1_lin_response` folder. The linear response formalism for
the calculation of the Hubbard U is based on the application
of a series of small perturbations to the first magnetic
atom. Automag isolates the first atom of the structure by
fictitiously changing its type. For this you need that

- the magnetic atomic type is listed first in your input
POSCAR file
- the VASP pseudopotential library on your supercomputer,
for the fictitious atomic type written at line 30 of the
`submit.py` file contains the same POSCAR file as your
actual magnetic atomic type (make this substitution by hand).

Before launching `submit.py`, you can edit the location of
the input structure at line 16, the magnetic state to use
for spin initialization at line 19 and the VASP parameters
at lines 36-47. Then launch the script. It will submit a
single VASP calculation: the bare run.

When the calculation has ended on your supercomputer, put
its calculation folder at line 25 of `submit.py`, change
the mode to 'NSC' at line 22 and relaunch the script.
Finally, change the mode to 'SC' and launch it again.
The non-selfconsistent and the selfconsistent runs will
be submitted.

When calculations are done, from the output file in the
`fw_calcs` folder you can get the FireWorks ID of each
calculation and, using the command `lpad get_launchdir [ID]`
you can know its calculation folder. For each
non-selfconsistent and selfconsistent run, check from the
`OUTCAR` file that the magnetic moments have not changed
their orientation from the initialized values and copy the
values of the charge on the partially occupied electron
shell of the first atom in the file `plot_results.py` on
your local machine at lines 22 and 23, respectively. Double
check that they correspond to the perturbations listed at
line 19. Then lauch the `plot_results.py` script. A PNG
file will be produced with the linear interpolations of the
response functions for the non-selfconsistent and the
selfconsistent runs, with the values of their slopes
indicated on the graph. The value of the Hubbard U can be
obtained from U = 1/X - 1/X0, where X and X0 are the
selfconsistent and the non-selfconsistent slopes, respectively.

## Search for the most stable magnetic state

Enter the directory `2_coll` on your local machine and open
the file `run_enumlib.py`. You need to put the correct unit
cell composition at line 16, you need to specify which atoms
have non-zero magnetization at line 17 and to point to the
correct input file with the structure at line 23. At line 20,
you can write the maximum size of the supercell for generating
antiferromagnetic configurations. Then launch the script.
After its execution, you will see a new folder named `enumlib`
with the results of the run. The most important files in this
folder are the ones which start with the word 'vasp'. They all
contain your initial geometry in POSCAR format, but with the
magnetic atoms split into two groups: one to be assigned spin
up and the other to be assigned spin down.

At this point you can open the file `submit.py`, which is used
to save to the remote database a single point energy
calculation for each antiferromagnetic configuration in the
`enumlib` folder, as well as an energy calculation for the
ferromagnetic state. You can put the correct unit cell
composition at line 18, the magnetic atomic type at line 19,
the values used to initialize magnetic moments at lines 20-21
and the VASP parameters at lines 27-53. Then launch the script.

When all the single point energy calculations are done on
your supercomputer, copy the output files that you find in the
`fw_calcs` folder to the `2_coll` folder on your local machine
and open the script `plot_results.py`. Put the correct
filenames at line 16, the magnetic atomic type at line 17 and
run the script. It will print a warning on screen if some
calculations did not converge and it will produce the file
`stability.png`, which is a histogram of the final energies of
all magnetic configurations which have been calculated. If the
magnetic moments did not change orientation after the energy
calculation the corresponding bar is blue, otherwise it is red.

## Calculation of the critical temperature

Automag can calculate the critical temperature of the
magnetically ordered to paramagnetic phase transition from a
Monte Carlo simulation of the effective Hamiltonian which
describes the magnetic interaction. For this you need to have
the program VAMPIRE installed in your system.

Open the file `get_configurations.py` in the folder
`3_monte_carlo`, tune the input parameters at lines 20-27 and
run the script. For the most abundant geometry among those
generated by enumlib, it will save the magnetic configurations
and their corresponding energies in the output files
`configurations.txt` and `energies.txt`. Then open the script
`coupling_constants.py` and tune the input parameters at lines
17-21. It will read the output files of the previous script and,
for all nearest neighbors in the magnetic sublattice at a
distance less than the given cutoff, it will calculate the
coupling constants of the corresponding Heisemberg model from
a least squares interpolation of the arising system of
equations. In addition, it will produce an output file
`model.png` with a comparison of the Heisemberg model energies
with the DFT energies, which allows you to estimate the accuracy
of the model. The script `run_vampire.py` produces the input
file `vamp.ucf` for the software VAMPIRE, which is error-prone
and time-consuming to write by hand, based on the output from
the previous script. For this you need to report the
information printed on screen by `coupling_constants.py` at
lines 15-17 of `run_vampire.py` and to run the script.
A template of the other input files for the program VAMPIRE is
provided in the `vampire_input` folder.

Once the VAMPIRE calculation has ended, you can place the
`output` file written by VAMPIRE in the `3_monte_carlo` folder
and run the script `plot_results.py`. It will produce a file
`magnetization.png` which contains a plot of the mean
magnetization length with respect to temperature, from which
you can estimate the critical temperature of the magnetically
ordered to paramagnetic phase transition.

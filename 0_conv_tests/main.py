import os

from copy import copy
from common.calculators import VASPCalculator
from common.submission import Local


# initialization
flag = True
parameters_list = []
tags = ['ENCUT', 'KSPACING', 'SIGMA']

# collect input
while flag:
    print('Which VASP parameter do you want to check for convergence? (press a number or \'q\' for quit)')
    for i, tag in enumerate(tags):
        print(f'  {i + 1}. {tag}')

    inp = input()

    if inp.lower() == 'q':
        exit()

    index = int(inp) - 1
    if index < 0 or index > len(tags) - 1:
        raise ValueError(f'Value {tags} not understood.')

    if len(parameters_list) != 0 and tags[index] in parameters_list[0]:
        raise ValueError(f'{tags[index]} parameter already defined.')

    print(f'Insert the trial values for the {tags[index]} parameter separated by spaces:')

    inp = input()

    new_parameters_list = []
    if len(parameters_list) == 0:
        for value in inp.split():
            new_parameters_list.append({tags[index]: value})
    else:
        for parameters in parameters_list:
            for value in inp.split():
                new_parameters = copy(parameters)
                new_parameters[tags[index]] = value
                new_parameters_list.append(new_parameters)

    parameters_list = new_parameters_list

    print('Do you want to check another parameter simultaneously? (y/n)')

    inp = input()

    if inp.lower() != 'y':
        flag = False

# create root folder
rootfold = '-'.join([param.lower() for param in parameters_list[0].keys()])

if not os.path.isdir(rootfold):
    os.mkdir(rootfold)

calculator = VASPCalculator(os.path.join(os.getcwd(), 'Specific'))
submission = Local()

for parameters in parameters_list:
    calcfold = '-'.join(parameters.values())
    calculator.prepare(parameters, calcfold)
    submission.submit('mpirun -n 4 vasp_std', calcfold)

import os

from common.calculators import VASPCalculator

# initialize calculator
calculator = VASPCalculator(os.path.join(os.getcwd(), 'Specific'))

for value in range(250, 300, 10):
    parameters = {'encut': value}
    calculator.submit(parameters, os.path.join(os.getcwd(), 'encut', str(value)))

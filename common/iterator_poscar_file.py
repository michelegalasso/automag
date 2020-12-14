"""
@file:      iterator_poscar_file.py
@author:    Sergey Pozdnyakov
@brief:     Iterator for files of type gatheredPOSCARS.
"""


def iterator_poscar_file(file_name):
    with open(file_name, 'r') as f:
        while True:
            current_line = f.readline()
            if current_line == '':
                return
            lines = [current_line]
            for i in range(6):
                lines.append(f.readline())
            atom_nums = lines[-1].split()
            num_atoms = 0
            for i in range(len(atom_nums)):
                num_atoms += int(atom_nums[i])

            for i in range(num_atoms + 1):
                lines.append(f.readline())
            total_line = ""
            for line in lines:
                total_line += line
            yield total_line[0:-1]

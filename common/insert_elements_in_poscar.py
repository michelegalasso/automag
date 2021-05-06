"""
automag.common.insert_elements_in_poscar
========================================

Simple function which inserts element names in poscars from enumlib.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""


def insert_elements_in_poscar(poscar_string, composition, magnetic_atom):
    to_replace = '  '
    replace_with_1 = ''
    replace_with_2 = ''
    for element, abundance in composition.items():
        if element.name == magnetic_atom:
            to_replace += f'{int(abundance / 2)}   {int(abundance / 2)}   '
        else:
            to_replace += f'{int(abundance)}   '

        replace_with_1 += f'{element.name:>3s} '
        replace_with_2 += f'{int(abundance):3d} '

    return poscar_string.replace(to_replace, replace_with_1 + '\n' + replace_with_2)

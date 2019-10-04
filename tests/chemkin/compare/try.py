"""
compare thermo
"""

import os
import numpy as np
from chemkin_io import mechparser as mparser


# Get mechanism information
def _read_file(file_name):
    with open(file_name, encoding='utf8', errors='ignore') as file_obj:
        file_str = file_obj.read()
    return file_str


PATH = os.path.dirname(os.path.realpath(__file__))
MECH1_PATH = os.path.join(PATH, '../data/test')
MECH1_STR = _read_file(os.path.join(MECH1_PATH, 'm1.txt'))
MECH1_CSV_STR = _read_file(os.path.join(MECH1_PATH, 'm1.csv'))
MECH2_PATH = os.path.join(PATH, '../data/test')
MECH2_STR = _read_file(os.path.join(MECH2_PATH, 'm2.txt'))
MECH2_CSV_STR = _read_file(os.path.join(MECH2_PATH, 'm2.csv'))

INDEX = 'inchi'

# Temperatures and Pressures to run
T_REF = 1.0
TEMPS = np.array([500.0, 1000.0, 1500.0])
# TEMPS = [300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0,
#          1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0,
#          2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0]
PRESSURES = np.array([1.0, 4.0, 5.0])


def test__csvread():
    """ test 
    """

    # Build dictionaries containing:
    # thermo data strings for each species, k data strings for each reaction
    # Dictionaries indexed by the given mechanism names or InCHI string
    if INDEX == 'name':
        _, m2_thm_dct = mparser.compare.thermo.build_name_dcts(
            MECH1_STR, MECH2_STR)
        m1_rxn_dct, m2_rxn_dct = mparser.compare.rates.build_name_dcts(
            MECH1_STR, MECH2_STR)
    elif INDEX == 'inchi':
        _, m2_thm_dct = mparser.compare.thermo.build_inchi_dcts(
            MECH1_STR, MECH2_STR, MECH1_CSV_STR, MECH2_CSV_STR)
        m1_rxn_dct, m2_rxn_dct = mparser.compare.rates.build_inchi_dcts(
            MECH1_STR, MECH2_STR, MECH1_CSV_STR, MECH2_CSV_STR)

    print('\nMech1 reaction dct')
    for key, val in m1_rxn_dct.items():
        print(key)
        print(val)
    print('\n\nMech2 reaction dct')
    for key, val in m2_rxn_dct.items():
        print(key)
        print(val)

if __name__ == '__main__':
    test__csvread()

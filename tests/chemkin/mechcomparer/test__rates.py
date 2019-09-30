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


def test__compare_rates():
    """ test chemkin_io.mechparser.compare.rates
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

    # Calculate units to keep all comparisons straight...
    m1_units = mparser.mechanism.reaction_units(MECH1_STR)
    m2_units = mparser.mechanism.reaction_units(MECH2_STR)
    print('\n\nMech1 Units')
    print(m1_units)
    print('\n\nMech2 Units')
    print(m2_units)

    # Rate constant
    ktp_dct = mparser.compare.rates.calculate_reaction_rates(
        m1_rxn_dct, m2_rxn_dct,
        m1_units, m2_units,
        m2_thm_dct,
        T_REF, TEMPS, PRESSURES)

    # print dict
    print('\n\nktp dct')
    for rxn, mechs in ktp_dct.items():
        print(rxn)
        m1, m2 = mechs['m1'], mechs['m2']
        for (p1, k1), (p2, k2) in zip(m1.items(), m2.items()):
            print('Pressure: ', p1, p2)
            print('Mech1 ks: ', k1)
            print('Mech2 ks: ', k2)
            print(' ')


if __name__ == '__main__':
    test__compare_rates()

"""
compare thermo
"""

import os
from chemkin_io import mechparser as mparser


RC = 1.98720425864083e-3  # Gas Constant in kcal/mol.K


# Get mechanism information
def _read_file(file_name):
    with open(file_name, encoding='utf8', errors='ignore') as file_obj:
        file_str = file_obj.read()
    return file_str


PATH = os.path.dirname(os.path.realpath(__file__))

MECH1_PATH = os.path.join(PATH, '../mechparser/data/syngas')
MECH1_STR = _read_file(os.path.join(
    MECH1_PATH, 'mechanism_green.txt'))
MECH1_CSV_STR = _read_file(os.path.join(
    MECH1_PATH, 'smiles_green.csv'))

MECH2_PATH = os.path.join(PATH, '../mechparser/data/syngas')
MECH2_STR = _read_file(os.path.join(
    MECH2_PATH, 'mechanism_elte_2016_ijck.txt'))
MECH2_CSV_STR = _read_file(os.path.join(
    MECH2_PATH, 'smiles_elte_alt.csv'))

INDEX = 'inchi'

# Temperatures to run
# TEMPS = [500.0, 1000.0, 1500.0]
TEMPS = [300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0,
         1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0,
         2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0]


def test__compare_thermo():
    """ test chemkin_io.mechparser.compare.thermo
    """

    # Build dictionaries containing the thermo data strings for each species
    # Dictionaries indexed by the given mechanism names or InCHI string
    if INDEX == 'name':
        m1_thm_dct, m2_thm_dct = mparser.compare.thermo.build_name_dcts(
            MECH1_STR, MECH2_STR)
    elif INDEX == 'inchi':
        m1_thm_dct, m2_thm_dct = mparser.compare.thermo.build_inchi_dcts(
            MECH1_STR, MECH2_STR, MECH1_CSV_STR, MECH2_CSV_STR)

    print('\nMech1 thermo dct')
    for key, val in m1_thm_dct.items():
        print(key)
        print(val)
    print('\n\nMech2 thermo dct')
    for key, val in m2_thm_dct.items():
        print(key)
        print(val)

    # Calculate Enthalpy, Entropy, Gibbs, & Heat Capacity for each species
    # Mech 1 and 2 thermo data collated for each species
    thm_vals_dct = mparser.compare.thermo.calculate_mech_thermo(
        m1_thm_dct, m2_thm_dct, TEMPS)

    print('\n\n\nCombined thermo vals')
    for idx in thm_vals_dct:

        # Obtain the name of each species to print for ID purposes
        if INDEX == 'name':
            m1_name = idx
            m2_name = idx
        elif INDEX == 'inchi':
            print('\n\nInCHI: ', idx)
            m1_name, m2_name = mparser.compare.therm.get_mech_name_for_species(
                MECH1_CSV_STR, MECH2_CSV_STR, idx)
        print('M1 Name: ', m1_name)
        print('M2 Name: ', m2_name)

        # Print all of the thermo quantities
        m1_vals = thm_vals_dct[idx]['m1']
        m2_vals = thm_vals_dct[idx]['m2']
        print('\nM1 Enthalpy', m1_vals[0])
        print('M2 Enthalpy', m2_vals[0])
        print('M1 Entropy', m1_vals[1])
        print('M2 Entropy', m2_vals[1])
        print('M1 Gibbs', m1_vals[2])
        print('M2 Gibbs', m2_vals[2])
        print('M1 Heat Capacity', m1_vals[3])
        print('M2 Heat Capacity', m2_vals[3])


if __name__ == '__main_':
    test__compare_thermo()

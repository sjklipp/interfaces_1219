"""
compare thermo
"""

import os
from chemkin_io import mechparser as mparser


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
    MECH2_PATH, 'mechanism_elte_2016_ijck_alt.txt'))
MECH2_CSV_STR = _read_file(os.path.join(
    MECH2_PATH, 'smiles_elte_alt.csv'))

INDEX = 'inchi'

# Temperatures and Pressures to run
T_REF = 1.0
TEMPS = [500.0, 1000.0, 1500.0]
# TEMPS = [300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0,
#          1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0,
#          2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0]
PRESSURES = [1.0, 5.0, 10.0]


def test__compare_rates():
    """ test chemkin_io.mechparser.compare.rates
    """

    # Build dictionaries containing:
    # thermo data strings for each species, k data strings for each reaction
    # Dictionaries indexed by the given mechanism names or InCHI string
    if INDEX == 'name':
        m1_thm_dct, m2_thm_dct = mparser.compare.thermo.build_name_dcts(
            MECH1_STR, MECH2_STR)
        m1_rxn_dct, m2_rxn_dct = mparser.compare.rates.build_name_dcts(
            MECH1_STR, MECH2_STR)
    elif INDEX == 'inchi':
        m1_thm_dct, m2_thm_dct = mparser.compare.thermo.build_inchi_dcts(
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
    m1_units = mparser.compare.rate.reaction_units(MECH1_STR)
    m2_units = mparser.compare.rate.reaction_units(MECH2_STR)

    # Rate constant
    ktp_dct = mparser.compare.rates.calculate_reaction_rates(
        m1_rxn_dct, m2_rxn_dct,
        m1_units, m2_units,
        m2_thm_dct,
        T_REF, TEMPS, PRESSURES)

    # 


if __name__ == '__main__':
    test__compare_rates()

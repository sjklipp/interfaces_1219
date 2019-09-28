"""
compare thermo
"""

import os
import chemkin_io


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
TEMPS = [300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0,
         1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0,
         2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0]


# Functions to build dictionaries
def build_name_dcts(mech1_str, mech2_str):
    """ builds the thermo dictionaries indexed by names
    """

    mech1_thermo_block = chemkin_io.mechparser.util.clean_up_whitespace(
        chemkin_io.mechparser.mechanism.thermo_block(mech1_str))
    mech1_thermo_dct = chemkin_io.mechparser.thermo.dct_name_idx(
        mech1_thermo_block)

    mech2_thermo_block = chemkin_io.mechparser.util.clean_up_whitespace(
        chemkin_io.mechparser.mechanism.thermo_block(mech2_str))
    mech2_thermo_dct = chemkin_io.mechparser.thermo.dct_name_idx(
        mech2_thermo_block)

    return mech1_thermo_dct, mech2_thermo_dct


def build_inchi_dcts(mech1_str, mech2_str,
                     mech1_csv_str, mech2_csv_str):
    """ builds new thermo dictionaries indexed by inchis
    """

    mech1_thermo_block = chemkin_io.mechparser.util.clean_up_whitespace(
        chemkin_io.mechparser.mechanism.thermo_block(mech1_str))
    mech1_name_dct = chemkin_io.mechparser.mechanism.species_name_inchi_dct(
        mech1_csv_str)
    mech1_thermo_dct = chemkin_io.mechparser.thermo.dct_inchi_idx(
        mech1_thermo_block, mech1_name_dct)

    mech2_thermo_block = chemkin_io.mechparser.util.clean_up_whitespace(
        chemkin_io.mechparser.mechanism.thermo_block(mech2_str))
    mech2_name_dct = chemkin_io.mechparser.mechanism.species_name_inchi_dct(
        mech2_csv_str)
    mech2_thermo_dct = chemkin_io.mechparser.thermo.dct_inchi_idx(
        mech2_thermo_block, mech2_name_dct)

    return mech1_thermo_dct, mech2_thermo_dct


def calculate_mech_thermo(m1_thermo_dct, m2_thermo_dct, temps):
    """ Loop over the the Mech1 thermo entries
    """

    thermo_dct = {}

    # Calculate all thermo quanties with mech1 keys existing
    for m1_name, m1_thermo_dstr in m1_thermo_dct.items():

        # Calculate the thermo values for mech1
        m1_enthalpy, m1_entropy, m1_gibbs = [], [], []
        for temp in temps:
            m1_enthalpy.append(
                chemkin_io.mechparser.thermo.calculate_enthalpy(
                    m1_thermo_dstr, temp))
            m1_entropy.append(
                chemkin_io.mechparser.thermo.calculate_entropy(
                    m1_thermo_dstr, temp))
            m1_gibbs.append(
                chemkin_io.mechparser.thermo.calculate_gibbs(
                    m1_thermo_dstr, temp))

        # Calculate the thermo values for mech2
        m2_enthalpy, m2_entropy, m2_gibbs = [], [], []
        if m1_name in m2_thermo_dct:
            m2_thermo_dstr = m2_thermo_dct[m1_name]
            for temp in temps:
                m2_enthalpy.append(
                    chemkin_io.mechparser.thermo.calculate_enthalpy(
                        m2_thermo_dstr, temp))
                m2_entropy.append(
                    chemkin_io.mechparser.thermo.calculate_entropy(
                        m2_thermo_dstr, temp))
                m2_gibbs.append(
                    chemkin_io.mechparser.thermo.calculate_gibbs(
                        m2_thermo_dstr, temp))
        else:
            m2_enthalpy, m2_entropy, m2_gibbs = None, None, None

        # Add entry to overal thermo dictionary
        thermo_dct[m1_name] = {
            'm1': [m1_enthalpy, m1_entropy, m1_gibbs],
            'm2': [m2_enthalpy, m2_entropy, m2_gibbs],
        }

    # Now add the entries where m2 exists, but m1 does not
    uni_m2_names = [name
                    for name in m2_thermo_dct.keys()
                    if name not in m1_thermo_dct]
    for name in uni_m2_names:

        # Set values for mech1
        m1_enthalpy, m1_entropy, m1_gibbs = None, None, None

        # Calculate values for mech2
        m2_enthalpy, m2_entropy, m2_gibbs = [], [], []
        for temp in temps:
            m2_enthalpy.append(
                chemkin_io.mechparser.thermo.calculate_enthalpy(
                    m2_thermo_dct[name], temp))
            m2_entropy.append(
                chemkin_io.mechparser.thermo.calculate_entropy(
                    m2_thermo_dct[name], temp))
            m2_gibbs.append(
                chemkin_io.mechparser.thermo.calculate_gibbs(
                    m2_thermo_dct[name], temp))

        # Add entry to overal thermo dictionary
        thermo_dct[name] = {
            'm1': [m1_enthalpy, m1_entropy, m1_gibbs],
            'm2': [m2_enthalpy, m2_entropy, m2_gibbs],
        }


if __name__ == '__main__':

    if INDEX == 'name':
        M1_THERMO_DCT, M2_THERMO_DCT = build_name_dcts(
            MECH1_STR, MECH2_STR)
    elif INDEX == 'inchi':
        M1_THERMO_DCT, M2_THERMO_DCT = build_inchi_dcts(
            MECH1_STR, MECH2_STR, MECH1_CSV_STR, MECH2_CSV_STR)

    print('\nMech1 thermo dct')
    for key, val in M1_THERMO_DCT.items():
        print(key)
        print(val)
    print('\n\nMech2 thermo dct')
    for key, val in M2_THERMO_DCT.items():
        print(key)
        print(val)

    # calculate_mech_thermo(M1_THERMO_DCT, M2_THERMO_DCT, TEMPS)

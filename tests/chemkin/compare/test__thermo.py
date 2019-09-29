"""
compare thermo
"""

import os
import numpy as np
import chemkin_io


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
#TEMPS = [500.0, 1000.0, 1500.0]
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


def get_mech_name_for_species(mech1_csv_str, mech2_csv_str, ich):
    """ build dictionaries to get the name for a given InCHI string
    """

    mech1_inchi_dct = chemkin_io.mechparser.mechanism.species_inchi_name_dct(
        mech1_csv_str)
    mech2_inchi_dct = chemkin_io.mechparser.mechanism.species_inchi_name_dct(
        mech2_csv_str)

    if ich in mech1_inchi_dct:
        mech1_name = mech1_inchi_dct[ich]
    else:
        mech1_name = 'Not in Mechanism'
    if ich in mech2_inchi_dct:
        mech2_name = mech2_inchi_dct[ich]
    else:
        mech2_name = 'Not in Mechanism'

    return mech1_name, mech2_name


def calculate_mech_thermo(m1_thermo_dct, m2_thermo_dct, temps):
    """ Loop over the the Mech1 thermo entries
    """

    thermo_dct = {}

    # Calculate all thermo quanties with mech1 keys existing
    for m1_name, m1_thermo_dstr in m1_thermo_dct.items():

        # Calculate the thermo values for mech1
        m1_enthalpy, m1_entropy, m1_gibbs, m1_cp = [], [], [], []
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
            m1_cp.append(
                chemkin_io.mechparser.thermo.calculate_heat_capacity(
                    m1_thermo_dstr, temp))

        # Calculate the thermo values for mech2
        m2_enthalpy, m2_entropy, m2_gibbs, m2_cp = [], [], [], []
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
                m2_cp.append(
                    chemkin_io.mechparser.thermo.calculate_heat_capacity(
                        m2_thermo_dstr, temp))
        else:
            m2_enthalpy, m2_entropy, m2_gibbs, m2_cp = None, None, None, None

        # Add entry to overal thermo dictionary
        thermo_dct[m1_name] = {
            'm1': [m1_enthalpy, m1_entropy, m1_gibbs, m1_cp],
            'm2': [m2_enthalpy, m2_entropy, m2_gibbs, m2_cp],
        }

    # Now add the entries where m2 exists, but m1 does not
    uni_m2_names = [name
                    for name in m2_thermo_dct.keys()
                    if name not in m1_thermo_dct]
    for name in uni_m2_names:

        # Set values for mech1
        m1_enthalpy, m1_entropy, m1_gibbs, m1_cp = None, None, None, None

        # Calculate values for mech2
        m2_enthalpy, m2_entropy, m2_gibbs, m2_cp = [], [], [], []
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
            m2_cp.append(
                chemkin_io.mechparser.thermo.calculate_heat_capacity(
                    m2_thermo_dct[name], temp))

        # Add entry to overal thermo dictionary
        thermo_dct[name] = {
            'm1': [m1_enthalpy, m1_entropy, m1_gibbs, m1_cp],
            'm2': [m2_enthalpy, m2_entropy, m2_gibbs, m2_cp],
        }

    return thermo_dct


if __name__ == '__main__':

    if INDEX == 'name':
        M1_THM_DCT, M2_THM_DCT = build_name_dcts(
            MECH1_STR, MECH2_STR)
    elif INDEX == 'inchi':
        M1_THM_DCT, M2_THM_DCT = build_inchi_dcts(
            MECH1_STR, MECH2_STR, MECH1_CSV_STR, MECH2_CSV_STR)

    print('\nMech1 thermo dct')
    for key, val in M1_THM_DCT.items():
        print(key)
        print(val)
    print('\n\nMech2 thermo dct')
    for key, val in M2_THM_DCT.items():
        print(key)
        print(val)

    THM_VALS_DCT = calculate_mech_thermo(
        M1_THM_DCT, M2_THM_DCT, TEMPS)

    print('\n\n\nCombined thermo vals')
    for idx in THM_VALS_DCT:
        if INDEX == 'name':
            m1_name = idx
            m2_name = idx
        elif INDEX == 'inchi':
            print('\n\nInCHI: ', idx)
            m1_name, m2_name = get_mech_name_for_species(
                MECH1_CSV_STR, MECH2_CSV_STR, idx)
        print('M1 Name: ', m1_name)
        print('M2 Name: ', m2_name)
        m1_vals = THM_VALS_DCT[idx]['m1']
        m2_vals = THM_VALS_DCT[idx]['m2']
        print('\nM1 Enthalpy', m1_vals[0])
        print('M2 Enthalpy', m2_vals[0])
        print('M1 Entropy', m1_vals[1])
        print('M2 Entropy', m2_vals[1])
        print('M1 Gibbs', m1_vals[2])
        print('M2 Gibbs', m2_vals[2])
        print('M1 Heat Capacity', m1_vals[3])
        print('M2 Heat Capacity', m2_vals[3])

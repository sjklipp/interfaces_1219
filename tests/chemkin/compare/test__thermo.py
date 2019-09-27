"""
compare thermo
"""

import os
import chemkin_io


def _read_file(file_name):
    with open(file_name, encoding='utf8', errors='ignore') as file_obj:
        file_str = file_obj.read()
    return file_str


# Path to mechanism files
PATH = os.path.dirname(os.path.realpath(__file__))

# Get mechanism 1 information
MECH1_PATH = os.path.join(PATH, '../mechparser/data/syngas')
MECH1_STR = _read_file(os.path.join(MECH1_PATH, 'mechanism.txt'))

MECH1_THERMO_BLOCK = chemkin_io.mechparser.util.clean_up_whitespace(
    chemkin_io.mechparser.mechanism.thermo_block(MECH1_STR))
MECH1_THERMO_DCT = chemkin_io.mechparser.thermo.dct_name_idx(
        MECH1_THERMO_BLOCK)

# Get mechanism 2 information
MECH2_PATH = os.path.join(PATH, '../mechparser/data/syngas')
MECH2_STR = _read_file(os.path.join(MECH2_PATH, 'mechanism.txt'))

MECH2_THERMO_BLOCK = chemkin_io.mechparser.util.clean_up_whitespace(
    chemkin_io.mechparser.mechanism.thermo_block(MECH2_STR))
MECH2_THERMO_DCT = chemkin_io.mechparser.thermo.dct_name_idx(
        MECH2_THERMO_BLOCK)

# Temperatures
TEMPS = [300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0,
         1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0,
         2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0]


# Loop over the the Mech1 thermo entries
for m1_name, m1_thermo_dstr in MECH1_THERMO_DCT.items():

    # Loop over temperatures
    for temp in TEMPS:

        # Calculate the thermo values for mech1
        m1_enthalpy = chemkin_io.mechparser.thermo.calculate_enthalpy(
            m1_thermo_dstr, temp)
        m1_entropy = chemkin_io.mechparser.thermo.calculate_entropy(
            m1_thermo_dstr, temp)
        m1_gibbs = chemkin_io.mechparser.thermo.calculate_gibbs(
            m1_thermo_dstr, temp)

        # Calculate the thermo values for mech2
        m2_thermo_dstr = MECH2_THERMO_DCT[m1_name]
        m2_enthalpy = chemkin_io.mechparser.thermo.calculate_enthalpy(
            m2_thermo_dstr, temp)
        m2_entropy = chemkin_io.mechparser.thermo.calculate_entropy(
            m2_thermo_dstr, temp)
        m2_gibbs = chemkin_io.mechparser.thermo.calculate_gibbs(
            m2_thermo_dstr, temp)




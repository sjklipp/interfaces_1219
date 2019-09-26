""" test chemkin_io.mechparser.thermo
"""

from __future__ import unicode_literals
from builtins import open
import os
import numpy as np
import chemkin_io


def _read_file(file_name):
    with open(file_name, encoding='utf8', errors='ignore') as file_obj:
        file_str = file_obj.read()
    return file_str


PATH = os.path.dirname(os.path.realpath(__file__))
NATGAS_PATH = os.path.join(PATH, 'data/natgas')
NATGAS_MECH_STR = _read_file(os.path.join(NATGAS_PATH, 'mechanism.txt'))

NATGAS_THERMO_BLOCK = chemkin_io.mechparser.util.clean_up_whitespace(
    chemkin_io.mechparser.mechanism.thermo_block(NATGAS_MECH_STR))
NATGAS_BLOCK_STRS = chemkin_io.mechparser.thermo.data_strings(
    NATGAS_THERMO_BLOCK)
SPECIES_IDX = 100
SPECIES_POLYNOMIAL = NATGAS_BLOCK_STRS[SPECIES_IDX]

TEMP1 = 500.0
TEMP2 = 1000.0

# print('\n\nSpecies Polynomial from Mechanism File')
# print(SPECIES_POLYNOMIAL)


def test__species_name():
    """ test chemkin_io.mechparser.thermo.species_name
    """
    ref_name = 'S(1041)'
    name = chemkin_io.mechparser.thermo.species_name(
        SPECIES_POLYNOMIAL)
    assert name == ref_name


def test__temperatures():
    """ test chemkin_io.mechparser.thermo.temperatures
    """
    ref_temps = (100.0, 5000.0, 851.17)
    temps = chemkin_io.mechparser.thermo.temperatures(
        SPECIES_POLYNOMIAL)
    assert temps == ref_temps


def test__low_coefficients():
    """ test chemkin_io.mechparser.thermo.low_coefficients
    """
    ref_low_cfts = (1.38316741, 0.0607301618, -6.87644932e-05, 4.63273342e-08,
                    -1.30767374e-11, -27488.4597, 23.3435776)
    low_cfts = chemkin_io.mechparser.thermo.low_coefficients(
        SPECIES_POLYNOMIAL)
    assert low_cfts == ref_low_cfts


def test__high_coefficients():
    """ test chemkin_io.mechparser.thermo.high_coefficients
    """
    ref_high_cfts = (8.16250741, 0.0288724022, -1.26242468e-05, 2.35789979e-09,
                     -1.62798977e-13, -28642.576, -8.27062078)
    high_cfts = chemkin_io.mechparser.thermo.high_coefficients(
        SPECIES_POLYNOMIAL)
    assert high_cfts == ref_high_cfts


def test__data_block():
    """ test chemkin_io.mechparser.thermo.data_block
    """
    ref_block = ('S(1041)',
                 (100.0, 5000.0, 851.17),
                 (1.38316741, 0.0607301618, -6.87644932e-05, 4.63273342e-08,
                  -1.30767374e-11, -27488.4597, 23.3435776),
                 (8.16250741, 0.0288724022, -1.26242468e-05, 2.35789979e-09,
                  -1.62798977e-13, -28642.576, -8.27062078))
    block_str = chemkin_io.mechparser.thermo.data_block(
        NATGAS_THERMO_BLOCK)
    assert ref_block == block_str[SPECIES_IDX]


def test__data_strings():
    """ test chemkin_io.mechparser.thermo.data_strings
    """
    thm_strs = chemkin_io.mechparser.thermo.data_strings(
        NATGAS_THERMO_BLOCK)
    assert len(thm_strs) == 130


def test__temp_common_default():
    """ test chemkin_io.mechparser.thermo.temp_common_default
    """
    ref_temp_common = 1000.0
    temp_common = chemkin_io.mechparser.thermo.temp_common_default(
        NATGAS_THERMO_BLOCK)
    assert ref_temp_common == temp_common


def test__calculate_enthalpy():
    """ test chemkin_io.mechparser.thermo.calculate_enthalpy
    """
    ref_ht1 = -42.58312043165988
    ref_ht2 = -19.266388540056756
    ht1 = chemkin_io.mechparser.thermo.calculate_enthalpy(
        SPECIES_POLYNOMIAL, TEMP1)
    ht2 = chemkin_io.mechparser.thermo.calculate_enthalpy(
        SPECIES_POLYNOMIAL, TEMP2)
    assert np.isclose(ref_ht1, ht1)
    assert np.isclose(ref_ht2, ht2)


def test__calculate_entropy():
    """ test chemkin_io.mechparser.thermo.calculate_entropy
    """
    ref_st1 = 0.11016051269318868
    ref_st2 = 0.14192480698218363
    st1 = chemkin_io.mechparser.thermo.calculate_entropy(
        SPECIES_POLYNOMIAL, TEMP1)
    st2 = chemkin_io.mechparser.thermo.calculate_entropy(
        SPECIES_POLYNOMIAL, TEMP2)
    assert np.isclose(ref_st1, st1)
    assert np.isclose(ref_st2, st2)


def test__calculate_gibbs():
    """ test chemkin_io.mechparser.thermo.calculate_gibbs
    """
    ref_gt1 = -97.66337677825422
    ref_gt2 = -161.1911955222404
    gt1 = chemkin_io.mechparser.thermo.calculate_gibbs(
        SPECIES_POLYNOMIAL, TEMP1)
    gt2 = chemkin_io.mechparser.thermo.calculate_gibbs(
        SPECIES_POLYNOMIAL, TEMP2)
    assert np.isclose(ref_gt1, gt1)
    assert np.isclose(ref_gt2, gt2)


if __name__ == '__main__':
    test__species_name()
    test__temperatures()
    test__low_coefficients()
    test__high_coefficients()
    test__data_block()
    test__temp_common_default()
    test__calculate_enthalpy()
    test__calculate_entropy()
    test__calculate_gibbs()

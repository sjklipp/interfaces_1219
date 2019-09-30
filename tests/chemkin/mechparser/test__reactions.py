""" test chemkin_io.mechparser.mechanism
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
SYNGAS_PATH = os.path.join(PATH, '../data/syngas')
SYNGAS_MECH_STR = _read_file(os.path.join(SYNGAS_PATH, 'mechanism.txt'))

SYNGAS_REACTION_BLOCK = chemkin_io.mechparser.util.clean_up_whitespace(
    chemkin_io.mechparser.mechanism.reaction_block(SYNGAS_MECH_STR))
SYNGAS_REACTION_STRS = chemkin_io.mechparser.reaction.data_strings(
    SYNGAS_REACTION_BLOCK)
REACTION = SYNGAS_REACTION_STRS[20]

# TEST_PATH = os.path.join(PATH, '../data')
# TEST_MECH_STR = _read_file(os.path.join(TEST_PATH, 'test.txt'))
#
# TEST_REACTION_BLOCK = chemkin_io.mechparser.util.clean_up_whitespace(
#     chemkin_io.mechparser.mechanism.reaction_block(TEST_MECH_STR))
# TEST_REACTION_STRS = chemkin_io.mechparser.reaction.data_strings(
#     TEST_REACTION_BLOCK)
# TEST_REACTION = TEST_REACTION_STRS[0]

TROE_REACTION = SYNGAS_REACTION_STRS[0]
LINDEMANN_REACTION = SYNGAS_REACTION_STRS[2]
CHEBYSHEV_REACTION = SYNGAS_REACTION_STRS[12]
HIGHP_REACTION = 'CO(1)+O(5)(+M)=CO2(12)(+M)     1.880e+11 0.000     2.430'
PLOG_REACTION = """HOCO<=>CO+OH         6.300E+032    -5.960   32470.0
PLOG/      0.0010     1.550E-008     2.930      8768.0/
PLOG/      0.0030     1.770E+003     0.340     18076.0/
PLOG/      0.0296     2.020E+013    -1.870     22755.0/
PLOG/      0.0987     1.680E+018    -3.050     24323.0/
PLOG/      0.2961     2.500E+024    -4.630     27067.0/
PLOG/      0.9869     4.540E+026    -5.120     27572.0/"""
# print('\n\nRate Data for a Reaction')
# print(REACTION)

T_REF = 1.0
TEMPS = np.array([500.0, 1000.0, 1500.0, 2000.0])
PRESSURES = np.array([1, 5, 10])
PRESSURES2 = np.array([0.0100, 0.0700, 0.987])

# print('\nhigh p')
# print(HIGHP_REACTION)
# print('\nlindemann')
# print(LINDEMANN_REACTION)
# print('\ntroe')
# print(TROE_REACTION)
# print('\nchebyshev')
# print(CHEBYSHEV_REACTION)
print('\nplog')
print(PLOG_REACTION)


def test__reactant_names():
    """ test chemkin_io.mechparser.reaction.reactant_names
    """
    names = chemkin_io.mechparser.reaction.reactant_names(
        TEST_REACTION)
    print('\nreactants')
    print(names)


def test__product_names():
    """ test chemkin_io.mechparser.reaction.product_names
    """
    names = chemkin_io.mechparser.reaction.product_names(
        TEST_REACTION)
    print('\nproducts')
    print(names)


def test__high_p_parameters():
    """ test chemkin_io.mechparser.reaction.high_p_parameters
    """
    params = chemkin_io.mechparser.reaction.high_p_parameters(
        LINDEMANN_REACTION)
    print('\nhigh-pressure parameters')
    print(params)


def test__low_p_parameters():
    """ test chemkin_io.mechparser.reaction.low_p_parameters
    """
    params = chemkin_io.mechparser.reaction.low_p_parameters(
        LINDEMANN_REACTION)
    print('\nlow-pressure parameters')
    print(params)


def test__troe_parameters():
    """ test chemkin_io.mechparser.reaction.troe_parameters
    """
    params = chemkin_io.mechparser.reaction.troe_parameters(
        TROE_REACTION)
    print('\nTroe parameters')
    print(params)


def test__chebyshev_parameters():
    """ test chemkin_io.mechparser.reaction.chebyshev_parameters
    """
    params = chemkin_io.mechparser.reaction.chebyshev_parameters(
        CHEBYSHEV_REACTION)
    print('\nChebyshev parameters')
    print(params)


def test__plog_parameters():
    """ test chemkin_io.mechparser.reaction.plog_parameters
    """
    params = chemkin_io.mechparser.reaction.plog_parameters(
        PLOG_REACTION)
    print('\nPLog parameters')
    print(params)


def test__buffer_enhance_factors():
    """ test chemkin_io.mechparser.reaction.buffer_enhance_factors
    """
    fct_dct = chemkin_io.mechparser.reaction.buffer_enhance_factors(
        LINDEMANN_REACTION)
    print('\nLow Pressure Buffer Enhhancement Factors')
    print(fct_dct)


def test__units():
    """ test chemkin_io.mechparser.reaction.data_strings
    """
    units = chemkin_io.mechparser.reaction.units(
        SYNGAS_REACTION_BLOCK)
    assert units == ('kcal/mole', 'moles')


def test__high_p_rate_constants():
    """ test chemkin_io.mechparser.reaction.calculate_rate_constants
        for a reaction with only high-pressure params
    """
    units = chemkin_io.mechparser.reaction.units(
        SYNGAS_REACTION_BLOCK)
    ktp_dct = chemkin_io.mechparser.reaction.calculate_rate_constants(
        HIGHP_REACTION, T_REF, units, TEMPS, pressures=None)
    print('\nhigh-pressure rate_constants')
    for key, val in ktp_dct.items():
        print(key)
        print(val)


def test__lindemann_rate_constants():
    """ test chemkin_io.mechparser.reaction.calculate_rate_constants
        for a reaction with high-pressure and low-pressure params
    """
    units = chemkin_io.mechparser.reaction.units(
        SYNGAS_REACTION_BLOCK)
    ktp_dct = chemkin_io.mechparser.reaction.calculate_rate_constants(
        LINDEMANN_REACTION, T_REF, units, TEMPS, pressures=PRESSURES)
    print('\nLindemann rate_constants')
    for key, val in ktp_dct.items():
        print(key)
        print(val)


def test__troe_rate_constants():
    """ test chemkin_io.mechparser.reaction.calculate_rate_constants
        for a reaction with only high-pressure, low-pressure, and Troe params
    """
    units = chemkin_io.mechparser.reaction.units(
        SYNGAS_REACTION_BLOCK)
    ktp_dct = chemkin_io.mechparser.reaction.calculate_rate_constants(
        TROE_REACTION, T_REF, units, TEMPS, pressures=PRESSURES)
    print('\nTroe rate_constants')
    for key, val in ktp_dct.items():
        print(key)
        print(val)


def test__chebyshev_rate_constants():
    """ test chemkin_io.mechparser.reaction.calculate_rate_constants
        for a reaction with only high-pressure and Chebyshev params
    """
    units = chemkin_io.mechparser.reaction.units(
        SYNGAS_REACTION_BLOCK)
    ktp_dct = chemkin_io.mechparser.reaction.calculate_rate_constants(
        CHEBYSHEV_REACTION, T_REF, units, TEMPS, pressures=PRESSURES)
    print('\nChebyshev rate_constants')
    for key, val in ktp_dct.items():
        print(key)
        print(val)


def test__plog_rate_constants():
    """ test chemkin_io.mechparser.reaction.calculate_rate_constants
        for a reaction with only high-pressure and PLog params
    """
    units = ('cal/mole', 'moles')
    ktp_dct = chemkin_io.mechparser.reaction.calculate_rate_constants(
        PLOG_REACTION, T_REF, units, TEMPS, pressures=PRESSURES2)
    print('\nPLog rate_constants')
    for key, val in ktp_dct.items():
        print(key)
        print(val)


def test__reactant_and_product_names():
    """ test chemkin_io.mechparser.reaction.reactant_and_product_names
    """
    names = chemkin_io.mechparser.reaction.reactant_and_product_names(
        SYNGAS_REACTION_BLOCK)
    print('\nreactant and product names')
    for name in names:
        print(name)


def test__data_strings():
    """ test chemkin_io.mechparser.reaction.data_strings
    """

    rxn_strs = chemkin_io.mechparser.reaction.data_strings(
        TEST_REACTION_BLOCK)
    print('\ndata strings')
    for string in rxn_strs:
        print(string)
    # assert len(rxn_strs) == 1678


def test__dct_name_idx():
    """ test chemkin_io.mechparser.reaction.dct_name_idx
    """
    rxn_dct = chemkin_io.mechparser.reaction.dct_name_idx(
        SYNGAS_REACTION_BLOCK)
    for i, (key, val) in enumerate(rxn_dct.items()):
        print('\n')
        print(key)
        print(val)
        if i == 20:
            break


if __name__ == '__main__':
    # test__reactant_names()
    # test__product_names()
    # test__high_p_parameters()
    # test__low_p_parameters()
    # test__troe_parameters()
    # test__chebyshev_parameters()
    # test__plog_parameters()
    # test__buffer_enhance_factors()
    # test__units()
    # test__high_p_rate_constants()
    # test__lindemann_rate_constants()
    # test__troe_rate_constants()
    # test__chebyshev_rate_constants()
    test__plog_rate_constants()
    # test__reactant_and_product_names()
    # test__data_strings()
    # test__dct_name_idx()

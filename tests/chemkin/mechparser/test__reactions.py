""" test chemkin_io.mechparser.mechanism
"""

from __future__ import unicode_literals
from builtins import open
import os
import chemkin_io


def _read_file(file_name):
    with open(file_name, encoding='utf8', errors='ignore') as file_obj:
        file_str = file_obj.read()
    return file_str


PATH = os.path.dirname(os.path.realpath(__file__))
NATGAS_PATH = os.path.join(PATH, 'data/natgas')
HEPTANE_PATH = os.path.join(PATH, 'data/heptane')
NATGAS_MECH_STR = _read_file(os.path.join(NATGAS_PATH, 'mechanism.txt'))

NATGAS_REACTION_BLOCK = chemkin_io.mechparser.util.clean_up_whitespace(
    chemkin_io.mechparser.mechanism.reaction_block(NATGAS_MECH_STR))
NATGAS_REACTION_STRS = chemkin_io.mechparser.reaction.data_strings(
    NATGAS_REACTION_BLOCK)
REACTION = NATGAS_REACTION_STRS[20]

print('\n\nRate Data for a Reaction')
print(REACTION)


def test__reactant_names():
    """ test chemkin_io.mechparser.reaction.reactant_names
    """
    names = chemkin_io.mechparser.reaction.reactant_names(
        REACTION)
    print('\nreactants')
    print(names)


def test__product_names():
    """ test chemkin_io.mechparser.reaction.product_names
    """
    names = chemkin_io.mechparser.reaction.product_names(
        REACTION)
    print('\nproducts')
    print(names)


def test__high_p_parameters():
    """ test chemkin_io.mechparser.reaction.high_p_parameters
    """
    params = chemkin_io.mechparser.reaction.high_p_parameters(
        REACTION)
    print('\nhigh-pressure parameters')
    print(params)


def test__reactant_and_product_names():
    """ test chemkin_io.mechparser.reaction.reactant_and_product_names
    """
    names = chemkin_io.mechparser.reaction.reactant_and_product_names(
        NATGAS_REACTION_BLOCK)
    print('\nreactant and product names')
    for name in names:
        print(name)


def test__data_strings():
    """ test chemkin_io.mechparser.reaction.data_strings
    """

    rxn_strs = chemkin_io.mechparser.reaction.data_strings(
        NATGAS_REACTION_BLOCK)
    print('\ndata strings')
    for string in rxn_strs:
        print(string)
    assert len(rxn_strs) == 1678


if __name__ == '__main__':
    test__reactant_names()
    test__product_names()
    test__high_p_parameters()
    test__reactant_and_product_names()
    test__data_strings()

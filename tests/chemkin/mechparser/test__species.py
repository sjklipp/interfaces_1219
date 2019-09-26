""" test chemkin_io.mechparser.mechanism
"""

from __future__ import unicode_literals
from builtins import open
import os
import new_chemkin_io


def _read_file(file_name):
    with open(file_name, encoding='utf8', errors='ignore') as file_obj:
        file_str = file_obj.read()
    return file_str


PATH = os.path.dirname(os.path.realpath(__file__))
NATGAS_PATH = os.path.join(PATH, 'data/natgas')
NATGAS_MECH_STR = _read_file(os.path.join(NATGAS_PATH, 'mechanism.txt'))

NATGAS_SPECIES_BLOCK = new_chemkin_io.mechparser.util.clean_up_whitespace(
    new_chemkin_io.mechparser.mechanism.species_block(NATGAS_MECH_STR))

# print('\n\nSpecies Block from Mechanism File')
# print(NATGAS_SPECIES_BLOCK)


def test__species__names():
    """ test new_chemkin_io.mechparser.mechanism.species_block
        test new_chemkin_io.mechparser.species.names
    """

    spc_names = new_chemkin_io.mechparser.species.names(
        NATGAS_SPECIES_BLOCK)
    assert len(spc_names) == 130


if __name__ == '__main__':
    test__species__names()

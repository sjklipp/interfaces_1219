""" functions operating on the mechanism string
"""

from io import StringIO
import pandas
import autoparse.pattern as app
from automol.smiles import inchi
from chemkin_io.mechparser import util


def species_block(mech_str):
    """ species block
    """
    block_str = util.block(
        string=_clean_up(mech_str),
        start_pattern=app.one_of_these(['SPECIES', 'SPEC']),
        end_pattern='END'
    )
    return block_str


def reaction_block(mech_str):
    """ reaction block
    """
    block_str = util.block(
        string=_clean_up(mech_str),
        start_pattern=app.one_of_these(['REACTIONS', 'REAC']),
        end_pattern='END'
    )
    return block_str


def thermo_block(mech_str):
    """ thermo block
    """
    block_str = util.block(
        string=_clean_up(mech_str),
        start_pattern=app.one_of_these(['THERMO ALL', 'THERM ALL', 'THER ALL',
                                        'THERMO', 'THERM', 'THER']),
        end_pattern='END'
    )
    return block_str


def species_name_inchi_dct(csv_str):
    """ build a dictionary of name idx and inchi entry
    """
    csv_file = StringIO(csv_str)
    data = pandas.read_csv(csv_file)

    spc_dct = {}
    ichs = [inchi(smiles) for smiles in data.smiles]
    spc_dct = dict(zip(data.name, ichs))

    return spc_dct


def _clean_up(mech_str):
    mech_str = util.remove_line_comments(
        mech_str, delim_pattern=app.escape('!'))
    mech_str = util.clean_up_whitespace(mech_str)
    return mech_str

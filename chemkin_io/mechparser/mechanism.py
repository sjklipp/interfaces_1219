""" functions operating on the mechanism string
"""

from io import StringIO
import pandas
import autoparse.pattern as app
from automol.smiles import inchi as _inchi
from automol.inchi import smiles as _smiles
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


def reaction_units(mech_str):
    """ reaction units
    """
    units = util.reaction_units(
        string=_clean_up(mech_str),
        start_pattern=app.one_of_these(['REACTIONS', 'REAC']),
        units_pattern=app.one_or_more(
            app.one_of_these([app.LETTER, app.escape('/')])),
    )
    return units


def species_name_inchi_dct(csv_str):
    """ build a dictionary of name idx and inchi entry
    """
    data = _read_csv(csv_str)

    spc_dct = {}
    if hasattr(data, 'InChI'):
        spc_dct = dict(zip(data.name, data.InChI))
    elif hasattr(data, 'SMILES'):
        smiles = [_inchi(smiles) for smiles in data.SMILES]
        spc_dct = dict(zip(data.name, smiles))
    else:
        print('No "InChI" or "SMILES" column in csv file')

    return spc_dct


def species_name_smiles_dct(csv_str):
    """ build a dictionary of name idx and inchi entry
    """
    data = _read_csv(csv_str)

    spc_dct = {}
    if hasattr(data, 'SMILES'):
        spc_dct = dict(zip(data.name, data.SMILES))
    elif hasattr(data, 'InChI'):
        ichs = [_smiles(ich) for ich in data.InChI]
        spc_dct = dict(zip(data.name, ichs))
    else:
        print('No "SMILES" or "InChI" column in csv file')

    return spc_dct


def species_name_mult_dct(csv_str):
    """ build a dictionary of name idx and inchi entry
    """
    data = _read_csv(csv_str)

    spc_dct = {}
    if hasattr(data, 'mult'):
        spc_dct = dict(zip(data.name, data.mult))
    else:
        print('No "mult" column in csv file')

    return spc_dct


def species_name_charge_dct(csv_str):
    """ build a dictionary of name idx and inchi entry
    """
    data = _read_csv(csv_str)

    spc_dct = {}
    if hasattr(data, 'charge'):
        spc_dct = dict(zip(data.name, data.charge))
    else:
        spc_dct = dict(zip(data.name, [0 for name in data.name]))

    return spc_dct


def species_name_sens_dct(csv_str):
    """ build a dictionary of name idx and inchi entry
    """
    data = _read_csv(csv_str)

    spc_dct = {}
    if hasattr(data, 'sens'):
        spc_dct = dict(zip(data.name, data.sens))
    else:
        spc_dct = dict(zip(data.name, [0. for name in data.name]))

    return spc_dct


def species_inchi_name_dct(csv_str):
    """ build a dictionary of inchi idx and name entry
    """
    data = _read_csv(csv_str)

    spc_dct = {}
    if hasattr(data, 'InChI'):
        spc_dct = dict(zip(data.name, data.InChI))
    elif hasattr(data, 'SMILES'):
        ichs = [_inchi(smiles) for smiles in data.SMILES]
        spc_dct = dict(zip(data.name, ichs))
    else:
        raise ValueError

    return spc_dct


def _read_csv(csv_str):
    """ read the csv file:
    """
    csv_file = StringIO(csv_str)
    data = pandas.read_csv(csv_file, comment='!', quotechar="'")
    return data


def _clean_up(mech_str):
    mech_str = util.remove_line_comments(
        mech_str, delim_pattern=app.escape('!'))
    mech_str = util.clean_up_whitespace(mech_str)
    return mech_str

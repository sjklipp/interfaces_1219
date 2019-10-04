""" functions operating on the thermo block string
"""


import autoparse.pattern as app
import autoparse.find as apf
from chemkin_io.mechanalyzer.parser import util


RC = 1.98720425864083e-3  # in kcal/mol.K


# Functions which act on the entire thermo block of mechanism file #
def data_block(block_str):
    """ find all thermo data
    """
    thm_dstr_lst = data_strings(block_str)
    thm_dat_lst = tuple(zip(
        map(species_name, thm_dstr_lst),
        map(temperatures, thm_dstr_lst),
        map(low_coefficients, thm_dstr_lst),
        map(high_coefficients, thm_dstr_lst)))
    return thm_dat_lst


def data_strings(block_str):
    """ thermo strings
    """
    headline_pattern = (
        app.LINE_START + app.not_followed_by(app.one_of_these(
            [app.DIGIT, app.PLUS, app.escape('=')])) +
        app.one_or_more(app.NONNEWLINE) +
        app.escape('1') + app.LINE_END
    )
    thm_strs = util.headlined_sections(
        string=block_str.strip(),
        headline_pattern=headline_pattern,
    )
    return thm_strs


def dct_name_idx(block_str):
    """ build a dictionary indexes by the species' CHEMKIN mechanism name
    """
    thm_dstr_lst = data_strings(block_str)
    thm_name_dct = {}
    for string in thm_dstr_lst:
        name = species_name(string)
        thm_name_dct[name] = string

    return thm_name_dct


def dct_inchi_idx(block_str, name_inchi_dct):
    """ build a dictionary indexed by the species' InCHI string
    """
    thm_name_dct = dct_name_idx(block_str)

    thm_inchi_dct = {}
    for name, thm_dstr in thm_name_dct.items():
        thm_inchi_dct[name_inchi_dct[name]] = thm_dstr

    return thm_inchi_dct


def temp_common_default(block_str):
    """ temperature defaults from the thermo block
    """
    pattern = (app.STRING_START +
               app.UNSIGNED_FLOAT + app.LINESPACES +
               app.capturing(app.UNSIGNED_FLOAT) + app.LINESPACES +
               app.UNSIGNED_FLOAT)
    capture = apf.first_capture(pattern, block_str)
    assert capture
    tmp_com_def = float(capture)
    return tmp_com_def


# Functions which act on the nasa polynomials of a single species #
def species_name(thm_dstr):
    """ get the species name from a thermo data string
    """
    pattern = app.STRING_START + app.capturing(app.one_or_more(app.NONSPACE))
    spc = apf.first_capture(pattern, thm_dstr)
    return spc


def temperatures(thm_dstr):
    """ get the common temperature from a thermo data string
    """
    headline = apf.split_lines(thm_dstr)[0]
    pattern = (app.LINESPACES + app.capturing(app.UNSIGNED_FLOAT) +
               app.LINESPACES + app.capturing(app.UNSIGNED_FLOAT) +
               app.LINESPACES + app.capturing(app.UNSIGNED_FLOAT))
    captures = apf.first_capture(pattern, headline)
    assert captures
    tmps = tuple(map(float, captures))
    return tmps


def low_coefficients(thm_dstr):
    """ get the low temperature thermo coefficients
    """
    capture_lst = apf.all_captures(app.EXPONENTIAL_FLOAT, thm_dstr)
    assert len(capture_lst) in (14, 15)
    cfts = tuple(map(float, capture_lst[7:14]))
    return cfts


def high_coefficients(thm_dstr):
    """ get the high temperature thermo coefficients
    """
    capture_lst = apf.all_captures(app.EXPONENTIAL_FLOAT, thm_dstr)
    assert len(capture_lst) in (14, 15)
    cfts = tuple(map(float, capture_lst[:7]))
    return cfts


def _coefficients_for_specific_temperature(thm_dstr, temp):
    """ return the set of coefficients of the polynomial (low or high)
        that should be used for a given temperature
    """

    temps = temperatures(thm_dstr)
    if temps[0] < temp < temps[1]:
        cfts = low_coefficients(thm_dstr)
    elif temps[1] < temp < temps[2]:
        cfts = high_coefficients(thm_dstr)
    else:
        cfts = None
        # raise ValueError('Temperature outside range of NASA polynomial')

    return cfts

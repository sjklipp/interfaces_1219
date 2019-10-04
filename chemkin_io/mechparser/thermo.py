""" functions operating on the thermo block string
"""


import numpy as np
import autoparse.pattern as app
import autoparse.find as apf
from chemkin_io.mechparser import util

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


# functions which calculate quantiies using data from the thermo section #
def mech_thermo(nasa_dct, temps):
    """ Loop over the the Mech1 thermo entries
    """

    # Calculate all thermo quanties with mech1 keys existing
    mech_thermo_dct = {}
    for name, thermo_dstr in nasa_dct.items():

        # Calculate the thermo values for mech1
        enthalpy, heat_capacity, entropy, gibbs, = [], [], [], []
        for temp in temps:
            enthalpy.append(
                mechparser.thermo.calculate_enthalpy(
                    thermo_dstr, temp))
            heat_capacity.append(
                mechparser.thermo.calculate_heat_capacity(
                    thermo_dstr, temp))
            entropy.append(
                mechparser.thermo.calculate_entropy(
                    thermo_dstr, temp))
            gibbs.append(
                mechparser.thermo.calculate_gibbs(
                    thermo_dstr, temp))

        mech_thermo_dct[name] = [enthalpy, heat_capacity, entropy, gibbs]

    return mech_thermo_dct


def calculate_enthalpy(thm_dstr, temp):
    """ Calculate the Enthalpy [H(T)] of a species using the
        coefficients of its NASA polynomial
    """

    cfts = _coefficients_for_specific_temperature(thm_dstr, temp)

    if cfts is not None:
        enthalpy = (
            cfts[0] +
            ((cfts[1] * temp) / 2.0) +
            ((cfts[2] * temp**2) / 3.0) +
            ((cfts[3] * temp**3) / 4.0) +
            ((cfts[4] * temp**4) / 5.0) +
            (cfts[5] / temp)
        )
        enthalpy *= (RC * temp)
    else:
        enthalpy = 0.0

    return enthalpy


def calculate_entropy(thm_dstr, temp):
    """ Calculate the Entropy [S(T)] of a species using the
        coefficients of its NASA polynomial
    """
    cfts = _coefficients_for_specific_temperature(thm_dstr, temp)

    if cfts is not None:
        entropy = (
            (cfts[0] * np.log(temp)) +
            (cfts[1] * temp) +
            ((cfts[2] * temp**2) / 2.0) +
            ((cfts[3] * temp**3) / 3.0) +
            ((cfts[4] * temp**4) / 4.0) +
            (cfts[6])
        )
        entropy *= RC
    else:
        entropy = 0.0

    return entropy


def calculate_gibbs(thm_dstr, temp):
    """ Calculate the Gibbs Free Energy [H(T)] of a species using the
        coefficients of its NASA polynomial
    """

    enthalpy = calculate_enthalpy(thm_dstr, temp)
    entropy = calculate_entropy(thm_dstr, temp)
    if enthalpy is not None and entropy is not None:
        gibbs = enthalpy - (entropy * temp)
    else:
        gibbs = None

    return gibbs


def calculate_heat_capacity(thm_dstr, temp):
    """ Calculate the Heat Capacity [Cp(T)] of a species using the
        coefficients of its NASA polynomial
    """
    cfts = _coefficients_for_specific_temperature(thm_dstr, temp)

    if cfts is not None:
        heat_capacity = (
            cfts[0] +
            (cfts[1] * temp) +
            (cfts[2] * temp**2) +
            (cfts[3] * temp**3) +
            (cfts[4] * temp**4)
        )
        heat_capacity *= RC
    else:
        heat_capacity = 0.0

    return heat_capacity

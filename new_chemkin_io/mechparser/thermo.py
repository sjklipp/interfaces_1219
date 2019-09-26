""" functions operating on the thermo block string
"""


import numpy as np
import autoparse.pattern as app
import autoparse.find as apf
from new_chemkin_io.mechparser import util

RCONST = 1.98720425864083e-3  # in kcal/mol.K


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
        raise ValueError('Temperature outside range of NASA polynomial')

    return cfts


# functions which calculate quantiies using data from the thermo section #
def calculate_enthalpy(thm_dstr, temp):
    """ Calculate the Enthalpy [H(T)] of a species using the
        coefficients of its NASA polynomial
    """

    cfts = _coefficients_for_specific_temperature(thm_dstr, temp)

    enthalpy = (
        cfts[0] +
        ((cfts[1] * temp) / 2.0) +
        ((cfts[2] * temp**2) / 3.0) +
        ((cfts[3] * temp**3) / 4.0) +
        ((cfts[4] * temp**4) / 5.0) +
        (cfts[5] / temp)
    )
    enthalpy *= (RCONST * temp)

    return enthalpy


def calculate_entropy(thm_dstr, temp):
    """ Calculate the Entropy [S(T)] of a species using the
        coefficients of its NASA polynomial
    """
    cfts = _coefficients_for_specific_temperature(thm_dstr, temp)

    entropy = (
        (cfts[0] * np.log(temp)) +
        (cfts[1] * temp) +
        ((cfts[2] * temp**2) / 2.0) +
        ((cfts[3] * temp**3) / 3.0) +
        ((cfts[4] * temp**4) / 4.0) +
        (cfts[6])
    )
    entropy *= RCONST

    return entropy


def calculate_gibbs(thm_dstr, temp):
    """ Calculate the Gibbs Free Energy [H(T)] of a species using the
        coefficients of its NASA polynomial
    """

    enthalpy = calculate_enthalpy(thm_dstr, temp)
    entropy = calculate_entropy(thm_dstr, temp)
    if enthalpy and entropy:
        gibbs = enthalpy - (entropy * temp)
    else:
        gibbs = None

    return gibbs
# def nasa_to_equilibrium_constant(reacs_coefs, prds_coefs, temp):
#     """ Calculate the equilibrium constant for a reaction using the
#         coefficients of the reactant and product NASA polynomials
#     """
#
#     # Calculation deltagibbs
#     gibbs_reacs = 0.0
#     for coefs in reacs_coefs:
#         gibbs_reacs += nasa_to_gibbs(coefs, temp)
#     gibbs_prds = 0.0
#     for coefs in prds_coefs:
#         gibbs_prds += nasa_to_gibbs(coefs, temp)
#
#     delta_gibbs = gibbs_prds - gibbs_reacs
#
#    # Calculate the Equilibrium Constant
#    equil_const = np.exp(-delta_gibbs / (R * temp))
#
#    return equil_const

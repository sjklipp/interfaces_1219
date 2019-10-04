""" functions operating on the reactions block string
"""


import itertools
from qcelemental import constants as qcc
import autoparse.pattern as app
import autoparse.find as apf
from autoparse import cast as ap_cast
from chemkin_io.mechanalyzer.parser import util


# Constants and Conversion factors
# NAVO = qcc.constants.avogadro_constant
NAVO = 6.0221409e+23
CAL2KCAL = qcc.conversion_factor('cal/mol', 'kcal/mol')
J2KCAL = qcc.conversion_factor('J/mol', 'kcal/mol')
KJ2KCAL = qcc.conversion_factor('kJ/mol', 'kcal/mol')
KEL2KCAL = qcc.conversion_factor('kelvin', 'kcal/mol')


# Various strings needed to parse the data sections of the Reaction block
CHEMKIN_ARROW = (app.maybe(app.escape('<')) + app.escape('=') +
                 app.maybe(app.escape('>')))
CHEMKIN_PLUS_EM = app.PLUS + 'M'
CHEMKIN_PAREN_PLUS_EM = app.escape('(') + app.PLUS + 'M' + app.escape(')')

SPECIES_NAME_PATTERN = (
    r'[^\s=+\-]' +
    app.zero_or_more(app.one_of_these(
        [app.LETTER, app.DIGIT, app.escape('(+)'), r'[#,()\-]',
         app.escape('['), app.escape(']')])) +
    app.zero_or_more(app.PLUS)
)
SPECIES_NAMES_PATTERN = app.series(
    app.padded(SPECIES_NAME_PATTERN), app.padded(app.PLUS))

REACTION_PATTERN = (SPECIES_NAMES_PATTERN + app.padded(CHEMKIN_ARROW) +
                    SPECIES_NAMES_PATTERN)
COEFF_PATTERN = (app.NUMBER + app.LINESPACES + app.NUMBER +
                 app.LINESPACES + app.NUMBER)

# Constants
NAVO = 6.02214076e23

# def all_rate_constants(block_str):
#     """ get the rate constants
#     """
#     rxn_strs = data_strings(block_str)
#     highp_k_lst = list(map(high_p_parameters, rxn_strs))
#
#     return highp_k_lst


def units(block_str):
    """ get the units for the rate parameters
    """
    units_str = block_str.strip().splitlines()[0]
    units_lst = units_str.split()
    if units_lst:
        ea_units = units_lst[0].lower()
        a_units = units_lst[1].lower()
    else:
        ea_units = 'cal/mole'
        a_units = 'moles'

    return ea_units, a_units


def data_strings(block_str):
    """ reaction strings
    """
    rxn_strs = util.headlined_sections(
        string=block_str.strip(),
        headline_pattern=CHEMKIN_ARROW,
    )
    return rxn_strs


def dct_name_idx(block_str):
    """ build a dictionary with the name dictionary
    """
    rxn_dstr_lst = data_strings(block_str)
    rxn_dct = {}
    for string in rxn_dstr_lst:
        rct_name = reactant_names(string)
        prd_name = product_names(string)
        key = (rct_name, prd_name)
        if key not in rxn_dct.keys():
            rxn_dct[key] = string
        else:
            rxn_dct[key] += '\n'+string

    return rxn_dct


def dct_inchi_idx(block_str, name_inchi_dct):
    """ build a dictionary indexed by rct and prd InCHI strings
    """
    rxn_names_dct = dct_name_idx(block_str)

    rxn_inchi_dct = {}
    for rxn_names, rxn_dstr in rxn_names_dct.items():

        [rct_names, prd_names] = rxn_names
        rct_ichs, prd_ichs = (), ()
        for rct in rct_names:
            rct_ichs += ((name_inchi_dct[rct]),)
        for prd in prd_names:
            prd_ichs += ((name_inchi_dct[prd]),)

        rxn_inchi_dct[(rct_ichs, prd_ichs)] = rxn_dstr

    return rxn_inchi_dct


# Functions which act on a single reaction #
def reactant_names(rxn_dstr):
    """ reactant species names
    """
    pattern = _first_line_pattern(
        rct_ptt=app.capturing(SPECIES_NAMES_PATTERN),
        prd_ptt=SPECIES_NAMES_PATTERN,
        coeff_ptt=COEFF_PATTERN
    )
    string = apf.first_capture(pattern, rxn_dstr)
    names = _split_reagent_string(string)
    return names


def product_names(rxn_dstr):
    """ product species names
    """
    pattern = _first_line_pattern(
        rct_ptt=SPECIES_NAMES_PATTERN,
        prd_ptt=app.capturing(SPECIES_NAMES_PATTERN),
        coeff_ptt=COEFF_PATTERN
    )
    string = apf.first_capture(pattern, rxn_dstr)
    names = _split_reagent_string(string)
    return names


def high_p_parameters(rxn_dstr):
    """ high-pressure parameters
    """
    pattern = _first_line_pattern(
        rct_ptt=SPECIES_NAMES_PATTERN,
        prd_ptt=SPECIES_NAMES_PATTERN,
        coeff_ptt=app.capturing(COEFF_PATTERN)
    )
    string_lst = apf.all_captures(pattern, rxn_dstr)
    vals = []
    for string in string_lst:
        vals += ap_cast(string.split())
    return vals


def low_p_parameters(rxn_dstr):
    """ low-pressure parameters
    """
    pattern = (
        'LOW' +
        app.zero_or_more(app.SPACE) + app.escape('/') +
        app.SPACES + app.capturing(app.NUMBER) +
        app.SPACES + app.capturing(app.NUMBER) +
        app.SPACES + app.capturing(app.NUMBER) +
        app.zero_or_more(app.SPACE) + app.escape('/')
    )
    params = apf.first_capture(pattern, rxn_dstr)
    if params is not None:
        params = [float(val) for val in params]
    return params


def troe_parameters(rxn_dstr):
    """ troe parameters
    """
    pattern = (
        'TROE' +
        app.zero_or_more(app.SPACE) + app.escape('/') +
        app.SPACES + app.capturing(app.NUMBER) +
        app.SPACES + app.capturing(app.NUMBER) +
        app.SPACES + app.capturing(app.NUMBER) +
        app.SPACES + app.maybe(app.capturing(app.NUMBER)) +
        app.zero_or_more(app.SPACE) + app.escape('/')
    )
    params = apf.first_capture(pattern, rxn_dstr)
    if params is not None:
        params = [float(val) for val in params]
    return params


def chebyshev_parameters(rxn_dstr):
    """ chebyshev parameters
    """
    temp_pattern = (
        'TCHEB' + app.zero_or_more(app.SPACE) + app.escape('/') +
        app.SPACES + app.capturing(app.FLOAT) +
        app.SPACES + app.capturing(app.FLOAT) +
        app.zero_or_more(app.SPACE) + app.escape('/')
    )
    pressure_pattern = (
        'PCHEB' + app.zero_or_more(app.SPACE) + app.escape('/') +
        app.SPACES + app.capturing(app.FLOAT) +
        app.SPACES + app.capturing(app.FLOAT) +
        app.zero_or_more(app.SPACE) + app.escape('/')
    )
    alpha_dimension_pattern = (
        'CHEB' + app.zero_or_more(app.SPACE) + app.escape('/') +
        app.SPACES + app.capturing(app.INTEGER) +
        app.SPACES + app.capturing(app.INTEGER) +
        app.zero_or_more(app.SPACE) + app.escape('/')
    )
    alpha_elements_pattern = (
        'CHEB' + app.zero_or_more(app.SPACE) + app.escape('/') +
        app.series(
            app.capturing(app.SPACES + app.capturing(app.EXPONENTIAL_FLOAT)),
            app.SPACES
        ) +
        app.zero_or_more(app.SPACE) + app.escape('/')
    )

    cheb_temps = apf.first_capture(temp_pattern, rxn_dstr)
    cheb_pressures = apf.first_capture(pressure_pattern, rxn_dstr)
    alpha_dims = apf.first_capture(alpha_dimension_pattern, rxn_dstr)
    alpha_elms = apf.all_captures(alpha_elements_pattern, rxn_dstr)
    if not alpha_elms:
        alpha_elms = None

    if all(vals is not None
           for vals in (cheb_temps, cheb_pressures, alpha_dims, alpha_elms)):
        cheb_temps = [float(val) for val in cheb_temps]
        cheb_pressures = [float(val) for val in cheb_pressures]
        alpha_dims = [int(val) for val in alpha_dims]
        alpha_elms = [list(map(float, row)) for row in alpha_elms]
        cheb_params = cheb_temps, cheb_pressures, alpha_dims, alpha_elms
    else:
        cheb_params = None

    return cheb_params


def plog_parameters(rxn_dstr):
    """ gets parameters associated with plog strings
    """
    pattern = (
        'PLOG' +
        app.zero_or_more(app.SPACE) + app.escape('/') +
        app.SPACES + app.capturing(app.NUMBER) +
        app.SPACES + app.capturing(app.NUMBER) +
        app.SPACES + app.capturing(app.NUMBER) +
        app.SPACES + app.capturing(app.NUMBER) +
        app.zero_or_more(app.SPACE) + app.escape('/')
    )
    params = apf.all_captures(pattern, rxn_dstr)
    if params:
        params = [list(map(float, row)) for row in params]
    else:
        params = None

    return params


def buffer_enhance_factors(rxn_dstr):
    """ get the factors of speed-up from bath gas
    """
    species_char = app.one_of_these([
        app.LETTER, app.DIGIT,
        app.escape('('), app.escape(')'),
        app.UNDERSCORE])
    species_name = app.one_or_more(species_char)

    # Get the line that could have the bath gas buffer enhancements
    bath_line_pattern = (
        _first_line_pattern(
            rct_ptt=SPECIES_NAMES_PATTERN,
            prd_ptt=SPECIES_NAMES_PATTERN,
            coeff_ptt=COEFF_PATTERN) + '\n' +
        app.capturing(app.LINE)
    )
    bath_string = apf.first_capture(bath_line_pattern, rxn_dstr)

    # Check if this line has bath gas factors or is for something else
    # If factors in string, get factors
    bad_strings = ('DUPLICATE', 'LOW', 'TROE', 'CHEB', 'PLOG')
    if (any(string in bath_string for string in bad_strings)
            and bath_string.strip() != ''):
        factors = None
    else:
        bath_string = '\n'.join(bath_string.strip().split())
        factor_pattern = (
            app.capturing(species_name) +
            app.escape('/') +
            app.capturing(app.NUMBER) +
            app.escape('/')
        )
        baths = apf.all_captures(factor_pattern, bath_string)
        factors = {}
        for bath in baths:
            factors[bath[0]] = float(bath[1])

    return factors


# helper functions #
def _first_line_pattern(rct_ptt, prd_ptt, coeff_ptt):
    return (app.STRING_START +
            rct_ptt + app.padded(CHEMKIN_ARROW) + prd_ptt +
            app.LINESPACES + coeff_ptt)


def _split_reagent_string(rgt_str):

    def _interpret_reagent_count(rgt_cnt_str):
        _pattern = (app.STRING_START + app.capturing(app.maybe(app.DIGIT)) +
                    app.capturing(app.one_or_more(app.NONSPACE)))
        cnt, rgt = apf.first_capture(_pattern, rgt_cnt_str)
        cnt = int(cnt) if cnt else 1
        rgts = (rgt,) * cnt
        return rgts

    rgt_str = apf.remove(app.LINESPACES, rgt_str)
    rgt_str = apf.remove(CHEMKIN_PAREN_PLUS_EM, rgt_str)
    rgt_str = apf.remove(CHEMKIN_PLUS_EM, rgt_str)
    pattern = app.PLUS + app.not_followed_by(app.PLUS)
    rgt_cnt_strs = apf.split(pattern, rgt_str)
    rgts = tuple(itertools.chain(*map(_interpret_reagent_count, rgt_cnt_strs)))

    return rgts

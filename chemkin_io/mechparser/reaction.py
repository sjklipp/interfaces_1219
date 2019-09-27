""" functions operating on the reactions block string
"""
import numpy as np
import itertools
import autoparse.pattern as app
import autoparse.find as apf
from autoparse import cast as ap_cast
from chemkin_io.mechparser import util
import ratefit

# Various strings needed to parse the data sections of the Reaction block
CHEMKIN_ARROW = (app.maybe(app.escape('<')) + app.escape('=') +
                 app.maybe(app.escape('>')))
CHEMKIN_PLUS_EM = app.PLUS + 'M'
CHEMKIN_PAREN_PLUS_EM = app.escape('(') + app.PLUS + 'M' + app.escape(')')

SPECIES_NAME_PATTERN = (
    r'[^\s=+\-]' +
    app.zero_or_more(app.one_of_these(
        [app.LETTER, app.DIGIT, app.escape('(+)'), r'[#,()\-]'])) +
    app.zero_or_more(app.PLUS)
)
SPECIES_NAMES_PATTERN = app.series(
    app.padded(SPECIES_NAME_PATTERN), app.padded(app.PLUS))

REACTION_PATTERN = (SPECIES_NAMES_PATTERN + app.padded(CHEMKIN_ARROW) +
                    SPECIES_NAMES_PATTERN)
COEFF_PATTERN = (app.NUMBER + app.LINESPACES + app.NUMBER +
                 app.LINESPACES + app.NUMBER)


# Functions which act on the entire thermo block of mechanism file #
#                               exclude_names=('OHV', 'CHV', 'CH(6)')):
def reactant_and_product_names(block_str,
                               exclude_names=()):
    """ reactants and products, by species_names
    """

    rxn_strs = data_strings(block_str)
    rct_names_lst = list(map(reactant_names, rxn_strs))
    prd_names_lst = list(map(product_names, rxn_strs))
    rxn_names_lst = tuple(filter(
        lambda x: not any(name in exclude_names for name in x[0] + x[1]),
        zip(rct_names_lst, prd_names_lst)))

    return rxn_names_lst


def all_rate_constants(block_str):
    """ get the rate constants
    """
    rxn_strs = data_strings(block_str)
    highp_k_lst = list(map(high_p_parameters, rxn_strs))

    return highp_k_lst


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

    cheb_temps = [float(val) for val in cheb_temps]
    cheb_pressures = [float(val) for val in cheb_pressures]
    alpha_dims = [int(val) for val in alpha_dims]
    alpha_elms = [list(map(float, row)) for row in alpha_elms]

    return cheb_temps, cheb_pressures, alpha_dims, alpha_elms


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
    params = [list(map(float, row)) for row in params]

    return params


def low_p_buffer_enhance_factors(rxn_dstr):
    """ get the factors of speed-up from bath gas
    """
    pattern = (
        _first_line_pattern(
            rct_ptt=SPECIES_NAMES_PATTERN,
            prd_ptt=SPECIES_NAMES_PATTERN,
            coeff_ptt=COEFF_PATTERN) +
        app.series(
            app.capturing(
                app.NONNEWLINE +
                app.escape('/') +
                app.capturing(app.NUMBER) +
                app.escape('/') +
                app.SPACES),
            app.SPACES)
    )
    all_factors = apf.first_capture(pattern, rxn_dstr)
    factor_lst = []
    for factor in all_factors.strip().split():
        tmp = factor.split('/')
        factor_lst.append([tmp[0], tmp[1]])
    return factor_lst


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


# calculator functions
def calculate_rate_constants(rxn_str, t_ref, temps, pressures=None):
    """ calculate the rate constant using the rxn_string
    """
    assert pressures is not None

    rate_constants = {}

    # Read the parameters from the reactions string
    highp_params = high_p_parameters(rxn_str)
    lowp_params = low_p_parameters(rxn_str)
    troe_params = troe_parameters(rxn_str)
    chebyshev_params = chebyshev_parameters(rxn_str)
    plog_params = plog_parameters(rxn_str)

    # Calculate high_pressure rates
    highp_ks = ratefit.fxns.arrhenius(highp_params, t_ref, temps)
    rate_constants['high'] = highp_ks

    # Calculate pressure-dependent rate constants based on discovered params
    # Either (1) Plog, (2) Chebyshev, (3) Lindemann, or (4) Troe
    pdep_dct = {}
    if plog_params:
        pdep_dct = _plog(plog_params, pressures, temps, t_ref)
    elif chebyshev_params:
        pdep_dct = _chebyshev(chebyshev_params, pressures, temps)
    elif lowp_params:
        lowp_ks = ratefit.fxns.arrhenius(lowp_params, t_ref, temps)
        if not troe_params:
            pdep_dct = ratefit.fxns.lindemann(
                highp_ks, lowp_ks, pressures, temps)
        else:
            pdep_dct = _troe(troe_params, highp_ks, lowp_ks, pressures, temps)
    if pdep_dct:
        for key, val in pdep_dct.items():
            rate_constants[key] = val

    return rate_constants


def _plog(plog_params, pressures, temps, t_ref):
    """ calc plog
    """
    plog_dct = {}
    for params in plog_params:
        plog_dct[params[0]] = params[1:]
    pdep_dct = ratefit.fxns.plog(plog_dct, pressures, temps, t_ref)
    return pdep_dct


def _chebyshev(chebyshev_params, pressures, temps):
    """ calc chebyshev
    """
    tmin = chebyshev_params[0][0]
    tmax = chebyshev_params[0][1]
    pmin = chebyshev_params[1][0]
    pmax = chebyshev_params[1][1]
    alpha = np.array(chebyshev_params[3])
    pdep_dct = ratefit.fxns.chebyshev(
        alpha, tmin, tmax, pmin, pmax, pressures, temps)
    return pdep_dct


def _troe(troe_params, highp_ks, lowp_ks, pressures, temps):
    """ calc troe
    """
    if len(troe_params) == 3:
        ts2 = None
    elif len(troe_params) == 4:
        ts2 = troe_params[3]
    pdep_dct = ratefit.fxns.troe(
        highp_ks, lowp_ks, pressures, temps,
        troe_params[0], troe_params[1], troe_params[2], ts2=ts2)
    return pdep_dct

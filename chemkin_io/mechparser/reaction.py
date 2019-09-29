""" functions operating on the reactions block string
"""


import itertools
import numpy as np
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

    print('split reagent')
    print(rgt_str)
    rgt_str = apf.remove(app.LINESPACES, rgt_str)
    rgt_str = apf.remove(CHEMKIN_PAREN_PLUS_EM, rgt_str)
    rgt_str = apf.remove(CHEMKIN_PLUS_EM, rgt_str)
    pattern = app.PLUS + app.not_followed_by(app.PLUS)
    rgt_cnt_strs = apf.split(pattern, rgt_str)
    rgts = tuple(itertools.chain(*map(_interpret_reagent_count, rgt_cnt_strs)))

    return rgts


# calculator functions
def calculate_rate_constants(rxn_str, t_ref, rxn_units, temps, pressures=None):
    """ calculate the rate constant using the rxn_string
    """
    rate_constants = {}

    # Read the parameters from the reactions string
    highp_params = high_p_parameters(rxn_str)
    lowp_params = low_p_parameters(rxn_str)
    troe_params = troe_parameters(rxn_str)
    chebyshev_params = chebyshev_parameters(rxn_str)
    plog_params = plog_parameters(rxn_str)
    print('\nlocated params')
    print(highp_params)
    print(lowp_params)
    print(troe_params)
    print(chebyshev_params)
    print(plog_params)

    # Calculate high_pressure rates
    highp_params = _update_params(highp_params, rxn_units)
    highp_ks = ratefit.fxns.arrhenius(highp_params, t_ref, temps)
    rate_constants['high'] = highp_ks

    # Calculate pressure-dependent rate constants based on discovered params
    # Either (1) Plog, (2) Chebyshev, (3) Lindemann, or (4) Troe
    # Update units if necessary
    if any(params is not None
           for params in (plog_params, chebyshev_params, lowp_params)):
        assert pressures is not None

    pdep_dct = {}
    if plog_params is not None:
        updated_plog_params = []
        for params in plog_params:
            updated_plog_params.append(_update_params(params, rxn_units))
        pdep_dct = _plog(updated_plog_params, pressures, temps, t_ref)

    elif chebyshev_params is not None:
        pdep_dct = _chebyshev(chebyshev_params, pressures, temps)

    elif lowp_params is not None:
        lowp_params = _update_params(lowp_params, rxn_units)
        lowp_ks = ratefit.fxns.arrhenius(lowp_params, t_ref, temps)
        if troe_params is not None:
            pdep_dct = _troe(troe_params, highp_ks, lowp_ks, pressures, temps)
        else:
            pdep_dct = ratefit.fxns.lindemann(
                highp_ks, lowp_ks, pressures, temps)

    if pdep_dct:
        for key, val in pdep_dct.items():
            rate_constants[key] = val

    return rate_constants


def _update_params(params, rxn_units):
    """ change the units if necessary
        only needed for highp, lowp, and plog
    """
    # Figure out converstion factors
    if rxn_units[0] == 'cal/mole':
        ea_conv_factor = 1000.0
    else:
        ea_conv_factor = 1.0

    if rxn_units[1] == 'molecules':
        a_conv_factor = NAVO
    else:
        a_conv_factor = 1.0

    # update units of params
    if params is not None:
        params[2] *= ea_conv_factor
        if len(params) == 6:
            params[5] *= ea_conv_factor

        params[0] *= a_conv_factor
        if len(params) == 6:
            params[3] *= a_conv_factor

    return params


def _plog(plog_params, pressures, temps, t_ref):
    """ calc plog
    """
    plog_dct = {}
    for params in plog_params:
        plog_dct[params[0]] = params[1:]
    pdep_dct = ratefit.fxns.plog(plog_dct, t_ref, pressures, temps)
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

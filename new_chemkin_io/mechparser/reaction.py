""" functions operating on the reactions block string
"""
import itertools
import autoparse.pattern as app
import autoparse.find as apf
from autoparse import cast as ap_cast
from new_chemkin_io.mechparser import util

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
    string = apf.first_capture(pattern, rxn_dstr)
    vals = ap_cast(string.split())
    return vals


def plog_parameters(rxn_dstr):
    """ gets parameters associated with plog strings
    """
    pattern = (
        'PLOG' +
        app.SPACES + app.escape('/') +
        app.SPACES + app.capturing(app.NUMBER) +
        app.SPACES + app.capturing(app.NUMBER) +
        app.SPACES + app.capturing(app.NUMBER) +
        app.SPACES + app.capturing(app.NUMBER) +
        app.SPACES + app.escape('/')
    )
    params = apf.all_captures(pattern, rxn_dstr)
    return params


def low_p_parameters(rxn_dstr):
    """ low-pressure parameters
    """
    pattern = (
        'LOW' +
        app.SPACES + app.escape('/') +
        app.SPACES + app.capturing(app.NUMBER) +
        app.SPACES + app.capturing(app.NUMBER) +
        app.SPACES + app.capturing(app.NUMBER) +
        app.SPACES + app.escape('/')
    )
    params = apf.first_capture(pattern, rxn_dstr)
    return params


def troe_parameters(rxn_dstr):
    """ troe parameters
    """
    pattern = (
        'TROE' +
        app.SPACES + app.escape('/') +
        app.SPACES + app.capturing(app.NUMBER) +
        app.SPACES + app.capturing(app.NUMBER) +
        app.SPACES + app.capturing(app.NUMBER) +
        app.SPACES + app.maybe(app.capturing(app.NUMBER)) +
        app.SPACES + app.escape('/')
    )
    params = apf.first_capture(pattern, rxn_dstr)
    return params


# def low_p_buffer_enhance_factors(rxn_dstr):
#     """ get the factors of speed-up from bath gas
#     """
#     pattern = (
#         _first_line_pattern(
#             rct_ptt=SPECIES_NAMES_PATTERN,
#             prd_ptt=SPECIES_NAMES_PATTERN,
#             coeff_ptt=COEFF_PATTERN) +
#         app.capturing(
#             app.series(
#                app.STUFF +
#                app.escape('/') +
#                app.FLOAT +
#                app.escape('/')
#             )
#         )
#     )
#     all_factors = apf.first_capture(pattern, rxn_dstr)
#     factor_lst = []
#     for factor in all_factors.strip().split():
#         tmp = factor.split('/')
#         factor_lst.append([tmp[0], tmp[1]])
#     return factor_lst


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

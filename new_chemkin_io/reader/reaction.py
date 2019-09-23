""" functions operating on the reactions block string
"""
import itertools
import autoparse.pattern as app
import autoparse.find as apf
from autoparse import cast as ap_cast
from chemkin_io import util

CHEMKIN_ARROW = (app.maybe(app.escape('<')) + app.escape('=') +
                 app.maybe(app.escape('>')))


def reactant_and_product_names(block_str,
                               exclude_names=('OHV', 'CHV', 'CH(6)')):
    """ reactants and products, by species_names
    """
    rxn_strs = data_strings(block_str)
    rct_names_lst = list(map(DataString.reactant_names, rxn_strs))
    prd_names_lst = list(map(DataString.product_names, rxn_strs))
    rxn_names_lst = tuple(filter(
        lambda x: not any(name in exclude_names for name in x[0] + x[1]),
        zip(rct_names_lst, prd_names_lst)))
    return rxn_names_lst


def data_strings(block_str):
    """ reaction strings
    """
    rxn_strs = util.headlined_sections(
        string=block_str.strip(),
        headline_pattern=CHEMKIN_ARROW,
    )
    return rxn_strs


class DataString():
    """ a reaction data string class """
    SPECIES_NAME_PATTERN = (
        r'[^\s=+\-]' + app.zero_or_more(app.one_of_these(
            [app.LETTER, app.DIGIT, app.escape('(+)'), r'[#,()\-]'])) +
        app.zero_or_more(app.PLUS)
    )
    SPECIES_NAMES_PATTERN = app.series(
        app.padded(SPECIES_NAME_PATTERN), app.padded(app.PLUS))

    REACTION_PATTERN = (SPECIES_NAMES_PATTERN + app.padded(CHEMKIN_ARROW) +
                        SPECIES_NAMES_PATTERN)
    COEFF_PATTERN = (app.NUMBER + app.LINESPACES + app.NUMBER +
                     app.LINESPACES + app.NUMBER)

    @staticmethod
    def _first_line_pattern(rct_ptt, prd_ptt, coeff_ptt):
        return (app.STRING_START +
                rct_ptt + app.padded(CHEMKIN_ARROW) + prd_ptt +
                app.LINESPACES + coeff_ptt)

    @classmethod
    def reactant_names(cls, rxn_str):
        """ reactant species names
        """
        pattern = cls._first_line_pattern(
            rct_ptt=app.capturing(cls.SPECIES_NAMES_PATTERN),
            prd_ptt=cls.SPECIES_NAMES_PATTERN,
            coeff_ptt=cls.COEFF_PATTERN
        )
        string = apf.first_capture(pattern, rxn_str)
        names = _split_reagent_string(string)
        return names

    @classmethod
    def product_names(cls, rxn_str):
        """ product species names
        """
        pattern = cls._first_line_pattern(
            rct_ptt=cls.SPECIES_NAMES_PATTERN,
            prd_ptt=app.capturing(cls.SPECIES_NAMES_PATTERN),
            coeff_ptt=cls.COEFF_PATTERN
        )
        string = apf.first_capture(pattern, rxn_str)
        names = _split_reagent_string(string)
        return names

    @classmethod
    def high_p_coefficients(cls, rxn_str):
        """ high-pressure coefficients
        """
        pattern = cls._first_line_pattern(
            rct_ptt=cls.SPECIES_NAMES_PATTERN,
            prd_ptt=cls.SPECIES_NAMES_PATTERN,
            coeff_ptt=app.capturing(cls.COEFF_PATTERN)
        )
        string = apf.first_capture(pattern, rxn_str)
        vals = ap_cast(string.split())
        return vals


# helpers
CHEMKIN_PLUS_EM = app.PLUS + 'M'
CHEMKIN_PAREN_PLUS_EM = app.escape('(') + app.PLUS + 'M' + app.escape(')')


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


# FROM old_chemkin_io #
# def reaction_data(mech_str):
#     """ find all reaction data
#     """
#     a_key, e_key = reaction_unit_names(mech_str)
#     rxn_dstr_lst = reaction_data_strings(mech_str)
# 
#     rxn_lst = list(map(reaction_data_reaction_name, rxn_dstr_lst))
#     arrh_lst = list(map(reaction_data_high_p_coeffs, rxn_dstr_lst))
#     arrh_lst = _convert_units(arrh_lst, a_key=a_key, e_key=e_key)
#     rxn_dat_lst = tuple(zip(rxn_lst, arrh_lst))
#     return rxn_dat_lst
# 
# 
# def _convert_units(arrh_lst, a_key, e_key):
#     a_unit_dct = dict(A_UNITS)
#     e_unit_dct = dict(E_UNITS)
#     a_lst, b_lst, e_lst = zip(*arrh_lst)
#     if a_key is not None:
#         assert a_key in a_unit_dct
#         a_lst = numpy.multiply(a_lst, a_unit_dct[a_key])
#     if e_key is not None:
#         assert e_key in e_unit_dct
#         e_lst = numpy.multiply(e_lst, e_unit_dct[e_key])
#     return list(zip(a_lst, b_lst, e_lst))
# 
# 
# def reaction_data_strings(mech_str):
#     """ find all reaction data strings
#     """
#     block_str = remove_blanks(reactions_block(mech_str))
#     headline_pattern = CHEMKIN_ARROW
#     rxn_dat_lst = _headlined_sections(headline_pattern, block_str)
#     return rxn_dat_lst
# 
# 
# def reaction_data_reaction_name(rxn_dstr):
#     """ get the reaction name from a reaction data string
#     """
#     headline = apf.split_lines(rxn_dstr)[0]
#     pattern = (app.LINESPACES + app.NUMBER + app.LINESPACES + app.NUMBER +
#                app.LINESPACES + app.NUMBER)
#     rxn = apf.remove(pattern, headline)
#     return rxn
# 
# 
# def split_reaction_name(rxn):
#     """ split a CHEMKIN reaction name into reactants and products
#     """
#     rct_str, prd_str = apf.split(CHEMKIN_ARROW, rxn)
#     rcts = _split_reagent_string(rct_str)
#     prds = _split_reagent_string(prd_str)
#     return rcts, prds
# def reaction_unit_names(mech_str):
#     """ units specified in the reaction block
#     """
#     block_str = remove_blanks(reactions_block(mech_str))
#     a_unit_names, _ = zip(*A_UNITS)
#     e_unit_names, _ = zip(*E_UNITS)
#     a_pattern = (app.STRING_START +
#                  app.maybe(app.one_of_these(e_unit_names) + app.LINESPACES) +
#                  app.capturing(app.one_of_these(a_unit_names)))
#     e_pattern = (app.STRING_START +
#                  app.maybe(app.one_of_these(a_unit_names) + app.LINESPACES) +
#                  app.capturing(app.one_of_these(e_unit_names)))
#     a_unit_name = apf.first_capture(a_pattern, block_str)
#     e_unit_name = apf.first_capture(e_pattern, block_str)
#     return a_unit_name, e_unit_name

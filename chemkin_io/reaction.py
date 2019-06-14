""" functions operating on the reactions block string
"""
import itertools
import autoparse.pattern as app
import autoparse.find as apf
from autoparse import cast as ap_cast
from chemkin_io import util

CHEMKIN_ARROW = (app.maybe(app.escape('<')) + app.escape('=') +
                 app.maybe(app.escape('>')))


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

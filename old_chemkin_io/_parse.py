""" CHEMKIN parsing

(could be cleaned up with new autoparse functionality)
"""
import itertools
import more_itertools as mit
import numpy
import autoparse.pattern as app
import autoparse.find as apf

CHEMKIN_ARROW = (app.maybe(app.escape('<')) + app.escape('=') +
                 app.maybe(app.escape('>')))
CHEMKIN_PLUS_EM = app.PLUS + 'M'
CHEMKIN_PAREN_PLUS_EM = app.escape('(') + app.PLUS + 'M' + app.escape(')')
A_UNITS = (
    ('MOLECULES', 6.02214076e23),
    ('MOLES', 1.)
)
E_UNITS = (
    ('KCAL/MOLE', 1.e3),
    ('CAL/MOLE', 1.),
    ('KJOULES/MOLE', 0.239006 * 1.e3),
    ('JOULES/MOLE', 0.239006),
    ('KELVINS', 0.001987191686485529 * 1.e3)
)


def species_names(mech_str):
    """ find all species
    """
    block_str = species_block(mech_str)
    spcs = apf.split_words(block_str)
    return spcs


def reaction_data(mech_str):
    """ find all reaction data
    """
    a_key, e_key = reaction_unit_names(mech_str)
    rxn_dstr_lst = reaction_data_strings(mech_str)

    rxn_lst = list(map(reaction_data_reaction_name, rxn_dstr_lst))
    arrh_lst = list(map(reaction_data_high_p_coeffs, rxn_dstr_lst))
    arrh_lst = _convert_units(arrh_lst, a_key=a_key, e_key=e_key)
    rxn_dat_lst = tuple(zip(rxn_lst, arrh_lst))
    return rxn_dat_lst


def _convert_units(arrh_lst, a_key, e_key):
    a_unit_dct = dict(A_UNITS)
    e_unit_dct = dict(E_UNITS)
    a_lst, b_lst, e_lst = zip(*arrh_lst)
    if a_key is not None:
        assert a_key in a_unit_dct
        a_lst = numpy.multiply(a_lst, a_unit_dct[a_key])
    if e_key is not None:
        assert e_key in e_unit_dct
        e_lst = numpy.multiply(e_lst, e_unit_dct[e_key])
    return list(zip(a_lst, b_lst, e_lst))


def reaction_data_strings(mech_str):
    """ find all reaction data strings
    """
    block_str = remove_blanks(reactions_block(mech_str))
    headline_pattern = CHEMKIN_ARROW
    rxn_dat_lst = _headlined_sections(headline_pattern, block_str)
    return rxn_dat_lst


def reaction_data_reaction_name(rxn_dstr):
    """ get the reaction name from a reaction data string
    """
    headline = apf.split_lines(rxn_dstr)[0]
    pattern = (app.LINESPACES + app.NUMBER + app.LINESPACES + app.NUMBER +
               app.LINESPACES + app.NUMBER)
    rxn = apf.remove(pattern, headline)
    return rxn


def split_reaction_name(rxn):
    """ split a CHEMKIN reaction name into reactants and products
    """
    rct_str, prd_str = apf.split(CHEMKIN_ARROW, rxn)
    rcts = _split_reagent_string(rct_str)
    prds = _split_reagent_string(prd_str)
    return rcts, prds


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


def reaction_data_high_p_coeffs(rxn_dstr):
    """ get the high-pressure Arrhenius coefficients from a reaction data string
    """
    headline = apf.split_lines(rxn_dstr)[0]
    pattern = (app.LINESPACES + app.capturing(app.NUMBER) + app.LINESPACES +
               app.capturing(app.NUMBER) + app.LINESPACES +
               app.capturing(app.NUMBER))
    captures = apf.first_capture(pattern, headline)
    assert captures
    cfts = tuple(map(float, captures))
    return cfts


def reaction_data_is_duplicate(rxn_dstr):
    """ is this a duplicate reaction?
    """
    pattern = 'DUP'
    return apf.has_match(pattern, rxn_dstr)


def thermo_data(mech_str):
    """ find all thermo data
    """
    thm_dstr_lst = thermo_data_strings(mech_str)
    thm_dat_lst = tuple(zip(
        map(thermo_data_species_name, thm_dstr_lst),
        map(thermo_data_lo_coefficients, thm_dstr_lst),
        map(thermo_data_hi_coefficients, thm_dstr_lst),
        map(thermo_data_temperatures, thm_dstr_lst)))
    return thm_dat_lst


def thermo_data_strings(mech_str):
    """ find all thermo data strings
    """
    block_str = remove_blanks(thermo_block(mech_str))
    start_pattern = app.LINE_START + app.not_followed_by(
        app.one_of_these([app.DIGIT, app.PLUS, app.escape('=')]))
    end_pattern = '1' + app.LINE_END
    headline_pattern = (start_pattern + app.one_or_more(app.NONNEWLINE) +
                        end_pattern)
    thm_dstr_lst = _headlined_sections(headline_pattern, block_str)
    assert all(len(apf.split_lines(thm_dstr)) == 4
               for thm_dstr in thm_dstr_lst)
    return thm_dstr_lst


def thermo_data_species_name(thm_dstr):
    """ get the species name from a thermo data string
    """
    pattern = app.STRING_START + app.capturing(app.one_or_more(app.NONSPACE))
    spc = apf.first_capture(pattern, thm_dstr)
    return spc


def thermo_data_temperatures(thm_dstr):
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


def thermo_data_lo_coefficients(thm_dstr):
    """ get the low temperature thermo coefficients
    """
    capture_lst = apf.all_captures(app.EXPONENTIAL_FLOAT, thm_dstr)
    assert len(capture_lst) in (14, 15)
    cfts = tuple(map(float, capture_lst[7:14]))
    return cfts


def thermo_data_hi_coefficients(thm_dstr):
    """ get the low temperature thermo coefficients
    """
    capture_lst = apf.all_captures(app.EXPONENTIAL_FLOAT, thm_dstr)
    assert len(capture_lst) in (14, 15)
    cfts = tuple(map(float, capture_lst[:7]))
    return cfts


def reaction_unit_names(mech_str):
    """ units specified in the reaction block
    """
    block_str = remove_blanks(reactions_block(mech_str))
    a_unit_names, _ = zip(*A_UNITS)
    e_unit_names, _ = zip(*E_UNITS)
    a_pattern = (app.STRING_START +
                 app.maybe(app.one_of_these(e_unit_names) + app.LINESPACES) +
                 app.capturing(app.one_of_these(a_unit_names)))
    e_pattern = (app.STRING_START +
                 app.maybe(app.one_of_these(a_unit_names) + app.LINESPACES) +
                 app.capturing(app.one_of_these(e_unit_names)))
    a_unit_name = apf.first_capture(a_pattern, block_str)
    e_unit_name = apf.first_capture(e_pattern, block_str)
    return a_unit_name, e_unit_name


def thermo_t_common_default(mech_str):
    """ temperature defaults from the thermo block
    """
    block_str = remove_blanks(thermo_block(mech_str))
    pattern = (app.STRING_START +
               app.UNSIGNED_FLOAT + app.LINESPACES +
               app.capturing(app.UNSIGNED_FLOAT) + app.LINESPACES +
               app.UNSIGNED_FLOAT)
    capture = apf.first_capture(pattern, block_str)
    assert capture
    tmp_com_def = float(capture)
    return tmp_com_def


def species_block(mech_str):
    """ find the species block
    """
    return _block(mech_str, block_keys=['SPECIES', 'SPEC'])


def reactions_block(mech_str):
    """ find the reactions block
    """
    return _block(mech_str, block_keys=['REACTIONS', 'REAC'])


def thermo_block(mech_str):
    """ find the thermodynamics block
    """
    return _block(mech_str, block_keys=['THERMO ALL', 'THERM ALL', 'THER ALL',
                                        'THERMO', 'THERM', 'THER'])


def remove_comments(mech_str):
    """ remove comments
    """
    pattern = app.escape('!') + app.zero_or_more(app.NONNEWLINE)
    clean_mech_str = apf.remove(pattern, mech_str)
    assert clean_mech_str is not None
    return clean_mech_str


def remove_blanks(mech_str):
    """ remove blank lines as well as leading and trailing blanks
    """
    blank_line = app.LINE_START + app.maybe(app.LINESPACES) + app.NEWLINE
    trailing_blanks = app.LINESPACES + app.LINE_END
    leading_blanks = app.LINE_START + app.LINESPACES
    pattern = app.one_of_these([blank_line, trailing_blanks, leading_blanks])
    return apf.remove(pattern, mech_str)


def _block(mech_str, block_keys):
    start_key = app.one_of_these(block_keys)
    contents = app.capturing(app.one_or_more(app.WILDCARD, greedy=False))
    pattern = start_key + contents + 'END'
    mech_str = remove_comments(mech_str)
    block = apf.first_capture(pattern, mech_str)
    return apf.strip_spaces(block) if block else None


def _headlined_sections(pattern, string):
    """ return sections with headlines matching a pattern
    """
    lines = string.splitlines()
    join_lines = '\n'.join
    pattern_matcher = apf.matcher(pattern)
    lines = mit.lstrip(lines, pred=lambda line: not pattern_matcher(line))
    sections = list(map(join_lines, mit.split_before(lines, pattern_matcher)))
    return sections

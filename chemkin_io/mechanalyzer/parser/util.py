""" utility functions
"""
import more_itertools as mit
import autoparse.pattern as app
import autoparse.find as apf


def clean_up_whitespace(string):
    """ remove leading spaces, trailing spaces, and empty lines from a string
    """
    empty_line = app.LINE_START + app.maybe(app.LINESPACES) + app.NEWLINE
    trailing_spaces = app.LINESPACES + app.LINE_END
    leading_spaces = app.LINE_START + app.LINESPACES
    pattern = app.one_of_these([empty_line, trailing_spaces, leading_spaces])
    return apf.remove(pattern, string)


def remove_line_comments(string, delim_pattern):
    """ remove line comments marked by a delimiter pattern
    """
    pattern = delim_pattern + app.zero_or_more(app.NONNEWLINE)
    return apf.remove(pattern, string)


def headlined_sections(string, headline_pattern):
    """ return sections with headlines matching a pattern
    """
    lines = string.splitlines()
    join_lines = '\n'.join
    pattern_matcher = apf.matcher(headline_pattern)
    lines = mit.lstrip(lines, pred=lambda line: not pattern_matcher(line))
    sections = list(map(join_lines, mit.split_before(lines, pattern_matcher)))
    return sections


def block(string, start_pattern, end_pattern):
    """ return a block delimited by start and end patterns
    """
    contents_pattern = app.capturing(
        app.one_or_more(app.WILDCARD, greedy=False))
    pattern = start_pattern + contents_pattern + end_pattern
    contents = apf.first_capture(pattern, string)
    return contents


def reaction_units(string, start_pattern, units_pattern):
    """ return a block delimited by start and end patterns
    """
    rxn_line_pattern = start_pattern + app.capturing(app.LINE_FILL)
    units_string = apf.first_capture(rxn_line_pattern, string)
    units_lst = apf.all_captures(units_pattern, units_string)

    ckin_ea_units = ['CAL/MOLE', 'KCAL/MOLE',
                     'JOULES/MOLE', 'KJOULES/MOLE',
                     'KELVINS']
    ckin_a_units = ['MOLES', 'MOLECULES']

    if units_lst:
        if any(unit in ckin_ea_units for unit in units_lst):
            for unit in ckin_ea_units:
                if unit in units_lst:
                    ea_unit = unit.lower()
        else:
            ea_unit = 'cal/mole'
        if any(unit in ckin_a_units for unit in units_lst):
            for unit in ckin_a_units:
                if unit in units_lst:
                    a_unit = unit.lower()
        else:
            a_unit = 'moles'
        units = (ea_unit, a_unit)
    else:
        units = ('cal/mole', 'moles')

    return units

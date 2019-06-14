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

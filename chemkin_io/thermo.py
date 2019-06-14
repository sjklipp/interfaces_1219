""" functions operating on the thermo block string
"""
import autoparse.pattern as app
from chemkin_io import util


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

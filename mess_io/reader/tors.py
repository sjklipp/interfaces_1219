"""
Reads the frequencies and ZPEs for torsional modes
from MESS output

Right now, we assume only a single species is in the output
"""

import autoparse.pattern as app
import autoparse.find as apf


def read_freqs(output_string):
    """ Reads the frequencies for the rotors
    """

    pattern = (app.escape('analytic  frequency at minimum[1/cm] =') +
               app.one_or_more(app.SPACE) +
               app.capturing(app.FLOAT))

    freqs = [float(val)
             for val in apf.all_captures(pattern, output_string)]

    return freqs


def read_zpes(output_string):
    """ Reads the ZPEs for the rotors
    """

    pattern = (app.escape('ground energy [kcal/mol]') +
               app.one_or_more(app.SPACE) +
               '=' +
               app.one_or_more(app.SPACE) +
               app.capturing(app.FLOAT))

    zpes = [float(val)
            for val in apf.all_captures(pattern, output_string)]

    return zpes

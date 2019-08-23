"""
Reads the outoput of a MESSPF calculation for the
frequencies and ZPEs for torsional modes

Right now, we assume only a single species is in the output
"""

import autoparse.pattern as app
import autoparse.find as apf


def freqs(output_string):
    """ Reads the frequencies for each hindered rotors from MESS output
        :param str output_string: string of lines for MESS output file
        :return freqs: List containing the frequencies for each rotor
        :rtype: list: float
    """

    # Pattern for the frequency of a rotor
    pattern = (app.escape('analytic  frequency at minimum[1/cm] =') +
               app.one_or_more(app.SPACE) +
               app.capturing(app.FLOAT))

    # Obtain each frequency from the output string
    tors_freqs = [float(val)
             for val in apf.all_captures(pattern, output_string)]

    return tors_freqs


def zpves(output_string):
    """ Reads the zero-point vibrational energies (ZPVES) for
        each hindered rotors from MESS output
        :param str output_string: string of lines for MESS output file
        :return freqs: List containing the ZPVES for each rotor
        :rtype: list: float
    """

    # Pattern for the ZPVE of a rotor
    pattern = (app.escape('ground energy [kcal/mol]') +
               app.one_or_more(app.SPACE) +
               '=' +
               app.one_or_more(app.SPACE) +
               app.capturing(app.FLOAT))

    # Obtain each ZPVE from the output string
    tors_zpes = [float(val)
            for val in apf.all_captures(pattern, output_string)]

    return tors_zpes

"""
Utility functions
"""

def indent(string, nspaces):
    """ Indents a multiline string. Required for Python 2,
        import textwrap textwrap.indent works for Python 3
    """
    pad = nspaces * ' '
    indented_string = ''.join(pad+line for line in string.splitlines(True))

    return indented_string


def elec_levels_format(elec_levels):
    """ Helps format the elec levels
    """

    # Get the number of elec levles
    nlevels = len(elec_levels)

    # Build elec levels string
    elec_levels_string = ''
    for i, level in enumerate(elec_levels):
        elec_levels_string += '  '.join(map(str, level))
        if (i+1) != len(elec_levels):
            elec_levels_string += '\n'

    # Indent the lines
    elec_levels_string = indent(elec_levels_string, 4)

    return nlevels, elec_levels_string


def geom_format(geom):
    """ Helps format the geometry
    """

    # Get the number of atoms
    natoms = len(geom)

    # Build geom string
    geom_string = ''
    for (asymb, xyz) in geom:
        geom_string += '{:s}   {:s}  {:s}  {:s}\n'.format(asymb, *map(repr, xyz))

    # Remove final newline character
    geom_string = geom_string.rstrip()

    # Indent the lines
    geom_string = indent(geom_string, 4)

    return natoms, geom_string


def freqs_format(freqs):
    """ Helps format the frequencies
    """

    # Get the number of freqs
    nfreqs = len(freqs)

    # Build freqs string
    freq_string = ''
    for i, freq in enumerate(freqs):
        if ((i+1) % 5) == 0 and (i+1) != len(freqs):
            freq_string += '{0:<8.2f}\n'.format(freq)
        else:
            freq_string += '{0:<8.2f}'.format(freq)

    # Indent the lines
    freq_string = indent(freq_string, 4)

    return nfreqs, freq_string


def format_rotor_key_defs(hind_rot_keyword_vals):
    """ formats Group, Axis, Symmetry keywords
    """

    # Build string containing the values of each keyword
    hind_rot_keyword_string = ''
    for vals in hind_rot_keyword_vals:
        hind_rot_keyword_string += '{0:<4d}'.format(vals)

    return hind_rot_keyword_string


def format_rotor_potential(potential):
    """ Helps format hindered rotor potential sections
    """

    # Get the number of the terms in the potential
    npotential = len(potential)

    # Build potentials string
    potential_string = ''
    for i, energy in enumerate(potential):
        if ((i+1) % 6) == 0 and (i+1) != npotential:
            potential_string += '{0:<8.2f}\n'.format(energy)
        else:
            potential_string += '{0:<8.2f}'.format(energy)

    # Indent the lines
    potential_string = indent(potential_string, 4)

    return npotential, potential_string

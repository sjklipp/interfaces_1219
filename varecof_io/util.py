"""
 Utility functions for formatting
"""

from qcelemental import constants as qcc
from qcelemental import periodictable as ptab
from automol import geom

# Conversion factors
BOHR2ANG = qcc.conversion_factor('bohr', 'angstrom')


def determine_struct_type(geo):
    """ determines the linear string
    """

    # Use automol to determine the type of structure
    if geom.is_atom(geo):
        struct_type = 'Monoatomic'
    else:
        if geom.is_linear(geo):
            struct_type = 'Linear'
        else:
            struct_type = 'Nonlinear'

    return struct_type


def format_coords(geo):
    """ format the coords section
    """

    # Get the number of atoms
    natoms = len(geo)

    # Get the geometry information
    symbols = geom.symbols(geo)
    coordinates = geom.coordinates(geo)
    masses = [int(ptab.to_mass(symbol)) for symbol in symbols]

    # Build a string with the formatted coordinates string
    if geom.is_atom(geo):
        geo_str = '{0:<4s}{1:<6d}'.format(symbols[0], masses[0])
    else:
        geo_str = ''
        for symbol, mass, coords in zip(symbols, masses, coordinates):
            coords = [coord * BOHR2ANG for coord in coords]
            coords_str = '{0:>14.8f}{1:>14.8f}{2:>14.8f}'.format(
                coords[0], coords[1], coords[2])
            geo_str += '{0:<4s}{1:<6d}{2}\n'.format(
                symbol, mass, coords_str)
        # Remove final newline character from the string
        geo_str = geo_str.rstrip()

    return natoms, geo_str

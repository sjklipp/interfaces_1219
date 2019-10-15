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

    # Remove dummy atoms
    geo = [coords for coords in geo
           if coords[0] != 'X']

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
    masses = [int(ptab.to_mass(symbol)) if symbol != 'X' else 0
              for symbol in symbols]

    # Build a string with the formatted coordinates string
    if geom.is_atom(geo):
        geo_str = '{0:<4s}{1:<6d}'.format(symbols[0], masses[0])
    else:
        geo_str = '{0} \n'.format(str(natoms))
        for symbol, mass, coords in zip(symbols, masses, coordinates):
            coords = [coord * BOHR2ANG for coord in coords]
            coords_str = '{0:>14.8f}{1:>14.8f}{2:>14.8f}'.format(
                coords[0], coords[1], coords[2])
            geo_str += '{0:<4s}{1:<6d}{2}\n'.format(
                symbol, mass, coords_str)
        # Remove final newline character from the string
        geo_str = geo_str.rstrip()

    return natoms, geo_str


def format_grids_string(grid, name, units):
    """ format the string using the grids for
        energy and angular momentum for tst.inp file
    """
    grid_str = '{0}_grid{1:>8d}{2:>9d}{3:>11.2f}{4:>7d}'.format(
        name, grid[0], grid[1], grid[2], grid[3])
    grid_str += '     {0:<8s}# {1} grid'.format(units, name)

    return grid_str


def format_faces_string(faces):
    """ format faces keyword
    """
    faces_str = ' '.join(faces)

    return faces_str


def format_values_string(coord, values, conv_factor=1.0):
    """ format the values string for the divsur.inp file
    """
    if values:
        values = ', '.join('{0:.3f}'.format(val * conv_factor)
                           for val in values)
        values_string = '{0} = ({1})'.format(coord, values)
    else:
        values_string = ''

    return values_string


def format_pivot_xyz_string(idx, npivot, xyzP, phi_dependence=False):
    """ format the pivot point xyz
    """

    assert npivot in (1, 2)

    atom_idx = idx
    if idx == 1:
        d_idx = 1
    else:
        d_idx = 2

    if npivot == 1:
        x_val = 'x{0} = {1:.3f}'.format(atom_idx, xyzP[0])
        y_val = '  y{0} = {1:.3f}'.format(atom_idx, xyzP[1])
        z_val = '  z{0} = {1:.3f}'.format(atom_idx, xyzP[2])
        pivot_xyz_string = (x_val + y_val + z_val)
    elif npivot > 1 and not phi_dependence:
        x_val1 = 'x{0} = {1:.3f} + d{2}*cos(t{3})'.format(
            atom_idx, xyzP[0], d_idx, atom_idx)
        y_val1 = '  y{0} = {1:.3f} + d{2}*sin(t{3})'.format(
            atom_idx, xyzP[1], d_idx, atom_idx)
        z_val1 = '  z{0} = 0.000'.format(
            atom_idx)
        x_val2 = 'x{0} = {1:.3f} - d{2}*cos(t{3})'.format(
            atom_idx+1, xyzP[0], d_idx, atom_idx)
        y_val2 = '  y{0} = {1:.3f} - d{2}*sin(t{3})'.format(
            atom_idx+1, xyzP[1], d_idx, atom_idx)
        z_val2 = '  z{0} = 0.000'.format(
            atom_idx+1)
        pivot_xyz_string = (x_val1 + y_val1 + z_val1 + '\n' +
                            x_val2 + y_val2 + z_val2)
    else:
        # Not sure if this implementation is any good
        x_val1 = 'x{0} = {1:.3f} + d{2}*sin(p{3})*cos(t{3})'.format(
            atom_idx, xyzP[0], d_idx, atom_idx)
        y_val1 = '  y{0} = {1:.3f} + d{2}*sin(p{3})*sin(t{3})'.format(
            atom_idx, xyzP[1], d_idx, atom_idx)
        z_val1 = '  z{0} = {1:.3f} + d{2}*cos(p{3})'.format(
            atom_idx, xyzP[2], d_idx, atom_idx)
        x_val2 = 'x{0} = {1:.3f} - d{2}*sin(p{3})*cos(t{3})'.format(
            atom_idx+1, xyzP[0], d_idx, atom_idx)
        y_val2 = '  y{0} = {1:.3f} - d{2}*sin(p{3})*sin(t{3})'.format(
            atom_idx+1, xyzP[1], d_idx, atom_idx)
        z_val2 = '  z{0} = {1:.3f} + d{2}*cos(p{3})'.format(
            atom_idx+1, xyzP[2], d_idx, atom_idx)
        pivot_xyz_string = (x_val1 + y_val1 + z_val1 + '\n' +
                            x_val2 + y_val2 + z_val2)

    return pivot_xyz_string

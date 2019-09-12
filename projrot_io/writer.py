"""
Functions to write the ProjRot input file
"""

import os
import numpy as np
from mako.template import Template
from qcelemental import constants as qcc
from qcelemental import periodictable as ptab


# Conversion factors
BOHR2ANG = qcc.conversion_factor('bohr', 'angstrom')

# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
TEMPLATE_PATH = os.path.dirname(os.path.realpath(__file__))


def rpht_input(geom, grad, hess,
               rotors_str='',
               coord_proj='cartesian'):
    """ Write the ProjRot input file
    """

    # Format the molecule info
    natoms = len(geom)
    geom_str = format_geom_str(geom)
    grad_str = format_grad_str(geom, grad)
    hess_str = format_hessian_str(hess)
    nrotors = rotors_str.count('pivotA')

    # Check what coordinate system to do projections in
    assert coord_proj in ('cartesian', 'internal')

    # Create a fill value dictionary
    rpht_keys = {
        'natoms': natoms,
        'nrotors': nrotors,
        'rotors_str': rotors_str,
        'geom_str': geom_str,
        'grad_str': grad_str,
        'hess_str': hess_str,
        'coord_proj': coord_proj
    }

    # Set template name and path for an atom
    template_file_name = 'rpht_input.mako'
    template_file_path = os.path.join(TEMPLATE_PATH, template_file_name)

    # Build a ProjRot input string
    rpht_string = Template(filename=template_file_path).render(**rpht_keys)

    return rpht_string


def rotors(axis, group):
    """ Write the sections that defines the rotors section
    """

    # Set up the keywords
    pivota = axis[0]
    pivotb = axis[1]
    atomsintopa = len(group)
    topaatoms = '  '.join([str(val) for val in group])

    # Build the rotors_str
    rotors_str = '{0:<32s}{1:<4d}\n'.format('pivotA', pivota)
    rotors_str += '{0:<32s}{1:<4d}\n'.format('pivotB', pivotb)
    rotors_str += '{0:<32s}{1:<4d}\n'.format('atomsintopA', atomsintopa)
    rotors_str += '{0:<32s}{1}'.format('topAatoms', topaatoms)

    return rotors_str


def format_geom_str(geo):
    """ Write the geometry section of the input file
        geometry in Angstroms
    """

    # Format the strings for the xyz coordinates
    geom_str = ''
    for i, (sym, coords) in enumerate(geo):
        anum = int(ptab.to_Z(sym))
        coords = [coord * BOHR2ANG for coord in coords]
        coords_str = '{0:>14.8f}{1:>14.8f}{2:>14.8f}'.format(
            coords[0], coords[1], coords[2])
        geom_str += '{0:2d}{1:4d}{2:4d}{3}\n'.format(
            i+1, anum, 0, coords_str)

    return geom_str


def format_grad_str(geom, grad):
    """ Write the gradient section of the input file
        grads in Hartrees/Bohr
    """

    atom_list = []
    for i, (sym, _) in enumerate(geom):
        atom_list.append(int(ptab.to_Z(sym)))

    # Format the strings for the xyz gradients
    full_grads_str = ''
    for i, grads in enumerate(grad):
        grads_str = '{0:>14.8f}{1:>14.8f}{2:>14.8f}'.format(
            grads[0], grads[1], grads[2])
        full_grads_str += '{0:2d}{1:4d}{2}\n'.format(
            i+1, atom_list[i], grads_str)

    return full_grads_str


def format_hessian_str(hess):
    """ Write the Hessian section of the input file
    """

    # Format the Hessian
    hess = np.array(hess)
    nrows, ncols = hess.shape

    if nrows % 5 == 0:
        nchunks = nrows // 5
    else:
        nchunks = (nrows // 5) + 1

    hess_str = '   ' + '    '.join([str(val) for val in range(1, 6)]) + '\n'
    cnt = 0
    while cnt+1 <= nchunks:
        for i in range(nrows):
            col_tracker = 1
            if i >= 5*cnt:
                hess_str += '{0}'.format(str(i+1))
                for j in range(5*cnt, ncols):
                    if i >= j:
                        if col_tracker <= 5:
                            hess_str += '  {0:5.8f}'.format(hess[i][j])
                            col_tracker += 1
                            if col_tracker == 6:
                                hess_str += '\n'
                        else:
                            continue
                    elif i < j and col_tracker != 6:
                        hess_str += '\n'
                        break
                    else:
                        break
            if i+1 == nrows and cnt+1 < nchunks-1:
                val_str = '     '.join(
                    [str(val)
                     for val in range(5*(cnt+1) + 1, 5*(cnt+1) + 6)])
                hess_str += '    ' + val_str + '\n'
            if i+1 == nrows and cnt+1 == nchunks-1:
                val_str = '     '.join(
                    [str(val)
                     for val in range(5*(cnt+1) + 1, nrows+1)])
                hess_str += '    ' + val_str + '\n'
        cnt += 1

    return hess_str

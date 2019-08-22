"""
Writes the global keyword section of a MESS input file
"""

import os
from mako.template import Template
from qcelemental import constants as qcc
from varecof_io import util


ANG2BOHR = qcc.conversion_factor('angstrom', 'bohr')

# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(SRC_PATH, 'templates')


def write_tst_input(nsamp_max, nsamp_min, flux_err, pes_size):
    """ Writes the tst.inp file for VaReCoF
        :param int nsamp_max: maximum number of samples
        :param int nsamp_min: minimum number of samples
        :return tst_inp_str: String for tst.inp file
        :rtype: string
    """

    # Create dictionary to fill template
    tst_keys = {
        'nsamp_max': nsamp_max,
        'nsamp_min': nsamp_min,
        'flux_err': flux_err,
        'pes_size': pes_size
    }

    # Set template name and path for the global keywords section
    template_file_name = 'tst.mako'
    template_file_path = os.path.join(TEMPLATE_PATH, template_file_name)

    # Build tst input string
    tst_str = Template(filename=template_file_path).render(**tst_keys)

    return tst_str


def write_divsur_input(distances):
    """ Writes the divsur.inp file for VaReCoF
        that contains info on the dividing surfaces.
        Right now we assume only center-of-mass seperations.
        :param float distances: List of temperatures (in Angstrom)
        :return divsur_inp_str: String for input file
        :rtype: string
    """

    # Format temperature and pressure lists
    distance_vals_str = ', '.join('{0:.3f}'.format(val * ANG2BOHR)
                                  for val in distances)
    distance_str = '(' + distance_vals_str + ')'

    # Create dictionary to fill template
    divsur_keys = {
        'distances': distance_str
    }

    # Set template name and path for the global keywords section
    template_file_name = 'divsur.mako'
    template_file_path = os.path.join(TEMPLATE_PATH, template_file_name)

    # Build divsur input string
    divsur_str = Template(filename=template_file_path).render(**divsur_keys)

    return divsur_str


def write_els_input(exe_path, base_name):
    """ Writes the electronic structure code input file for VaReCoF
        Currently code only runs with Molpro
        :rtype: string
    """

    # Create dictionary to fill template
    els_keys = {
        'exe_path': exe_path,
        'base_name': base_name
    }

    # Set template name and path for the global keywords section
    template_file_name = 'els.mako'
    template_file_path = os.path.join(TEMPLATE_PATH, template_file_name)

    # Build tst input string
    els_str = Template(filename=template_file_path).render(**els_keys)

    return els_str


def write_structure_input(geo1, geo2):
    """ Writes the structure input file for VaReCoF
        :rtype: string
    """

    # Determine linearity of molecule
    struct_type1 = util.determine_struct_type(geo1)
    struct_type2 = util.determine_struct_type(geo2)

    # Format the coordinates of the geoetry
    natoms1, coords1 = util.format_coords(geo1)
    natoms2, coords2 = util.format_coords(geo2)

    # Create dictionary to fill template
    struct_keys = {
        'struct_type1': struct_type1,
        'natoms1': natoms1,
        'coords1': coords1,
        'struct_type2': struct_type2,
        'natoms2': natoms2,
        'coords2': coords2,
    }

    # Set template name and path for the global keywords section
    template_file_name = 'struct.mako'
    template_file_path = os.path.join(TEMPLATE_PATH, template_file_name)

    # Build structure input string
    struct_str = Template(filename=template_file_path).render(**struct_keys)

    return struct_str


def write_tml_input(memory, basis, method, wfn, inf_sep_energy):
    """ writes the tml file used as the template for the electronic structure
        calculation
        currently, we assume the use of molpro
        in particular: method and wfn assume molpro input card structure
    """

    # convert the memory
    memory_mw = int(memory * (1024.0 / 8.0))

    # make the infinite seperation energy positive
    inf_sep_energy *= -1.0

    # Create dictionary to fill template
    tml_keys = {
        'memory': memory_mw,
        'basis': basis,
        'method': method,
        'wfn': wfn,
        'inf_sep_energy': inf_sep_energy
    }

    # Set template name and path for the tml input file string
    template_file_name = 'tml.mako'
    template_file_path = os.path.join(TEMPLATE_PATH, template_file_name)

    # Build structure input string
    tml_str = Template(filename=template_file_path).render(**tml_keys)

    return tml_str

"""
Writes MESS input for a molecule
"""

import os
from mako.template import Template
from mess_io.writer import util


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(SRC_PATH, 'templates')
SPECIES_PATH = os.path.join(TEMPLATE_PATH, 'species')
SPEC_INFO_PATH = os.path.join(SPECIES_PATH, 'info')


def write_core_rigidrotor(geom1, sym_factor, interp_emax=''):
    """ Writes a rigid-rotor core section.
    """

    # Format the geometry section
    natom1, geom1 = util.geom_format(geom1)

    # Create dictionary to fill template
    core_rigrot_keys = {
        'sym_factor': sym_factor,
        'natom1': natom1,
        'geom1': geom1,
        'interp_emax': interp_emax
    }

    # Set template name and path for a rigid-rotor core section
    template_file_name = 'core_rigidrotor.mako'
    template_file_path = os.path.join(SPEC_INFO_PATH, template_file_name)

    # Build core section string
    core_rigrot_str = Template(filename=template_file_path).render(**core_rigrot_keys)

    return core_rigrot_str


def write_core_multirotor(geom1, sym_factor, pot_surf, int_rot,
                          interp_emax=100, quant_lvl_emax=9):
    """ Writes a multi-rotor core section.
    """

    # Format the geometry section
    natom1, geom1 = util.geom_format(geom1)

    # Create dictionary to fill template
    core_multrot_keys = {
        'sym_factor': sym_factor,
        'natom1': natom1,
        'geom1': geom1,
        'pot_surf': pot_surf,
        'int_rot': int_rot,
        'interp_emax': interp_emax,
        'quant_lvl_emax': quant_lvl_emax
    }

    # Set template name and path for a multi-rotor core section
    template_file_name = 'core_multirotor.mako'
    template_file_path = os.path.join(SPEC_INFO_PATH, template_file_name)

    # Build core section string
    core_multrot_str = Template(filename=template_file_path).render(**core_multrot_keys)

    return core_multrot_str


def write_core_phasespace(geom1, geom2, sym_factor, stoich,
                          pot_prefactor=10, pot_power_exp=6):
    """ Writes a core section for phase space theory
    """

    # Format the geometry section of each fragment
    natom1, geom1 = util.geom_format(geom1)
    natom2, geom2 = util.geom_format(geom2)

    # Indent the geometry strings
    geom1 = util.indent(geom1, 2)
    geom2 = util.indent(geom2, 2)

    # Create dictionary to fill template
    core_phasespace_keys = {
        'sym_factor': sym_factor,
        'natom1': natom1,
        'geom1': geom1,
        'natom2': natom2,
        'geom2': geom2,
        'stoich': stoich,
        'pot_prefactor': pot_prefactor,
        'pot_power_exp': pot_power_exp
    }

    # Set template name and path for a phase space core section
    template_file_name = 'core_phasespace.mako'
    template_file_path = os.path.join(SPEC_INFO_PATH, template_file_name)

    # Build core section string
    core_phasespace_str = Template(filename=template_file_path).render(**core_phasespace_keys)

    return core_phasespace_str


def write_core_rotd(sym_factor, ne_file, stoich):
    """ Writes a core section which calls flux files from Rotd/VaReCoF
    """

    # Set values and template based on core type
    core_rotd_keys = {
        'sym_factor': sym_factor,
        'ne_file': ne_file,
        'stoich': stoich
    }

    # Set template name and path for a Rotd core section
    template_file_name = 'core_rotd.mako'
    template_file_path = os.path.join(SPEC_INFO_PATH, template_file_name)

    # Build core section string
    core_rotd_str = Template(filename=template_file_path).render(**core_rotd_keys)

    return core_rotd_str


def write_rotor_hindered(group, axis, symmetry, potential):
    """ Writes the section for a single hindered rotor.
    """

    # Format the sections
    rotor_hind_group = util.format_rotor_key_defs(group)
    rotor_hind_axis = util.format_rotor_key_defs(axis)
    rotor_hind_npotential, rotor_hind_potential = util.format_rotor_potential(potential)

    # Create dictionary to fill template
    rotor_hind_keys = {
        'group': rotor_hind_group,
        'axis': rotor_hind_axis,
        'symmetry': symmetry,
        'npotential': rotor_hind_npotential,
        'potential': rotor_hind_potential
    }

    # Set template name and path for a hindered rotor section
    template_file_name = 'rotor_hindered.mako'
    template_file_path = os.path.join(SPEC_INFO_PATH, template_file_name)

    # Build rotor section string
    rotor_hind_str = Template(filename=template_file_path).render(**rotor_hind_keys)

    return rotor_hind_str


def write_rotor_internal(group, axis, symmetry,
                         rotor_id='',
                         mass_exp_size=5, pot_exp_size=5,
                         hmin=13, hmax=101,
                         grid_size=100):
    """ Writes the section for a single internal rotor.
    """

    # Format the sections
    rotor_int_group = util.format_rotor_key_defs(group)
    rotor_int_axis = util.format_rotor_key_defs(axis)

    # Create dictionary to fill template
    rotor_int_keys = {
        'group': rotor_int_group,
        'axis': rotor_int_axis,
        'rotor_id': rotor_id,
        'symmetry': symmetry,
        'mass_exp_size': mass_exp_size,
        'pot_exp_size': pot_exp_size,
        'hmin': hmin,
        'hmax': hmax,
        'grid_size': grid_size
    }

    # Set template name and path for a internal rotor section
    template_file_name = 'rotor_internal.mako'
    template_file_path = os.path.join(SPEC_INFO_PATH, template_file_name)

    # Build rotor section string
    rotor_int_str = Template(filename=template_file_path).render(**rotor_int_keys)

    return rotor_int_str


def write_tunnel_eckart(imag_freq, well_depth1, well_depth2):
    """ Writes the tunneling section assuming an Eckart model
    """

    # Create dictionary to fill template
    tunnel_eckart_keys = {
        'imag_freq': imag_freq,
        'well_depth1': well_depth1,
        'well_depth2': well_depth2
    }

    # Set template name and path for an Eckart tunneling section
    template_file_name = 'tunnel_eckart.mako'
    template_file_path = os.path.join(SPEC_INFO_PATH, template_file_name)

    # Build tunnel section string
    tunnel_eckart_str = Template(filename=template_file_path).render(**tunnel_eckart_keys)

    return tunnel_eckart_str


def write_tunnel_sct(imag_freq, cutoff_energy, tunnel_file):
    """ Writes the tunneling section accounting for small curvature tunneling
    """

    # Create dictionary to fill template
    tunnel_sct_keys = {
        'imag_freq': imag_freq,
        'cutoff_energy': cutoff_energy,
        'tunnel_file': tunnel_file
    }

    # Set template name and path for an SCT tunneling section
    template_file_name = 'tunnel_sct.mako'
    template_file_path = os.path.join(SPEC_INFO_PATH, template_file_name)

    # Build tunnel section string
    tunnel_sct_str = Template(filename=template_file_path).render(**tunnel_sct_keys)

    return tunnel_sct_str

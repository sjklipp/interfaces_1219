"""
Writes MESS input for a molecule
"""

import os
from mako.template import Template
from mess_io.writer import util


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(SRC_PATH, 'templates')
SECTION_PATH = os.path.join(TEMPLATE_PATH, 'sections')
RXNCHAN_PATH = os.path.join(SECTION_PATH, 'reaction_channel')


def write_species(species_label, species_data):
    """ Writes a species section.
    """

    # Indent the string containing all of data for the well
    species_data = util.indent(species_data, 2)

    # Create dictionary to fill template
    species_keys = {
        'species_label': species_label,
        'species_data': species_data
    }

    # Set template name and path for a species
    template_file_name = 'species.mako'
    template_file_path = os.path.join(RXNCHAN_PATH, template_file_name)

    # Build species section string
    species_str = Template(filename=template_file_path).render(**species_keys)

    return species_str


def write_well(well_label, well_data):
    """ Writes a well section.
    """

    # Indent the string containing all of data for the well
    well_data = util.indent(well_data, 4)

    # Create dictionary to fill template
    well_keys = {
        'well_label': well_label,
        'well_data': well_data
    }

    # Set template name and path for a well
    template_file_name = 'well.mako'
    template_file_path = os.path.join(RXNCHAN_PATH, template_file_name)

    # Build well section string
    well_str = Template(filename=template_file_path).render(**well_keys)

    return well_str


def write_bimolecular(bimol_label,
                      species1_label, species1_data,
                      species2_label, species2_data,
                      ground_energy):
    """ Writes a Bimolecular section.
    """

    # Indent the string containing all of data for each species
    species1_data = util.indent(species1_data, 4)
    species2_data = util.indent(species2_data, 4)

    # Create dictionary to fill template
    bimol_keys = {
        'bimolec_label': bimol_label,
        'species1_label': species1_label,
        'species1_data': species1_data,
        'species2_label': species2_label,
        'species2_data': species2_data,
        'ground_energy': ground_energy
    }

    # Set template name and path for a bimolecular set 
    template_file_name = 'bimolecular.mako'
    template_file_path = os.path.join(RXNCHAN_PATH, template_file_name)

    # Build bimolecular section string
    bimolecular_str = Template(filename=template_file_path).render(**bimol_keys)

    return bimolecular_str


def write_ts_sadpt(ts_label, reac_label, prod_label, ts_data):
    """ Writes a TS section containing only a saddle point 
    """

    # Indent the string containing all of data for the saddle point
    ts_data = util.indent(ts_data, 2)

    # Create dictionary to fill template
    ts_sadpt_keys = {
        'ts_label': ts_label,
        'reac_label': reac_label,
        'prod_label': prod_label,
        'ts_data': ts_data
    }

    # Set template name and path for a TS with only a single saddle point
    template_file_name = 'ts_sadpt.mako'
    template_file_path = os.path.join(RXNCHAN_PATH, template_file_name)

    # Build saddle point string
    ts_sadpt_str = Template(filename=template_file_path).render(**ts_sadpt_keys)

    return ts_sadpt_str


def write_ts_irc(ts_label, reac_label, prod_label, irc_pt_strs):
    """ Writes a TS section containing IRC information
    """

    # Concatenate all of the IRC point strings
    ts_data = '\n'.join(irc_pt_strs)

    # Indent the TS IRC string
    ts_data = util.indent(ts_data, 4)

    # Create dictionary to fill template
    ts_irc_keys = {
        'ts_label': ts_label,
        'reac_label': reac_label,
        'prod_label': prod_label,
        'ts_data': ts_data
    }

    # Set template name and path for a TS with an IRC
    template_file_name = 'ts_irc.mako'
    template_file_path = os.path.join(RXNCHAN_PATH, template_file_name)

    # Build transition state with IRC string
    ts_irc_str = Template(filename=template_file_path).render(**ts_irc_keys)

    return ts_irc_str

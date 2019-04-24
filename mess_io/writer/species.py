"""
Writes MESS for an atom
"""

import os
from mako.template import Template
from writers import util


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(SRC_PATH, 'templates')
SPECIES_PATH = os.path.join(TEMPLATE_PATH, 'species')

def write_atom(atom_name, elec_levels):
    """ Writes the atom section of a MESS input file
    """

    # Build a formatted elec levels string
    nlevels, levels = util.elec_levels_format(elec_levels)

    # Create dictionary to fill template
    atom_keys = {
        'atom_name': atom_name,
        'nlevels': nlevels,
        'levels': levels
    }

    # Set template name and path for an atom
    template_file_name = 'atom.mako'
    template_file_path = os.path.join(SPECIES_PATH, template_file_name)

    # Build atom string
    atom_string = Template(filename=template_file_path).render(**atom_keys)

    return atom_string


def write_molecule(core, freqs, zero_energy, elec_levels,
                   hind_rot='', tunnel='',
                   anharm='', rovib_coups='', rot_dists=''):
    """ Write molecular info section
    """

    # Build a formatted frequencies and elec levels string
    nfreqs, freqs = util.freqs_format(freqs)
    nlevels, levels = util.elec_levels_format(elec_levels)

    # Format the rovib couplings and rotational distortions if needed
    if rovib_coups != '':
        rovib_coup_list = '  '.join(str(val) for val in rovib_coups)
    if rot_dists != '':
        rot_dists = util.format_rot_dist_consts(rot_dists)

    # Indent the anharm matrix string if needed
    if anharm != '':
        anharm = util.indent(anharm, 2)        

    # Create dictionary to fill template
    molec_keys = {
        'core': core,
        'nfreqs': nfreqs,
        'freqs': freqs,
        'zero_energy': zero_energy,
        'nlevels': nlevels,
        'levels': levels,
        'hind_rot': hind_rot,
        'anharm': anharm,
        'tunnel': tunnel
    }

    # Set template name and path for a molecule
    template_file_name = 'molecule.mako'
    template_file_path = os.path.join(SPECIES_PATH, template_file_name)

    # Build molecule string
    molecule_str = Template(filename=template_file_path).render(**molec_keys)

    return molecule_str

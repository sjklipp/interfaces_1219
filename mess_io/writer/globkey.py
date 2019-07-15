"""
Writes the global keyword section of a MESS input file
"""

import os
from mako.template import Template


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(SRC_PATH, 'templates')
SECTION_PATH = os.path.join(TEMPLATE_PATH, 'sections')

def write_global_reaction(temperatures, pressures):
    """ Writes the global keywords section of the MESS input file.
        :param float temperatures: List of temperatures (in K)
        :param float pressures: List of pressures (in atm)
        :return global_str: String for section
        :rtype: string
    """

    # Format temperature and pressure lists
    temperature_list = '  '.join(str(val) for val in temperatures)
    pressure_list = '  '.join(str(val) for val in pressures)

    # Create dictionary to fill template
    global_rxn_keys = {
        'temperatures': temperature_list,
        'pressures': pressure_list
    }

    # Set template name and path for the global keywords section
    template_file_name = 'global_reaction.mako'
    template_file_path = os.path.join(SECTION_PATH, template_file_name)

    # Build global section string
    global_rxn_str = Template(filename=template_file_path).render(**global_rxn_keys)

    return global_rxn_str


def write_global_pf(temp_step=100, ntemps=30, rel_temp_inc=0.001, atom_dist_min=0.6):
    """ Writes the global keywords section of the MESS input file.
        :param float temp_step: temperature step (in K)
        :param ntemps: number of temperature values on grid
        :return global_pf_str: String for section
        :rtype: string
    """

    # Create dictionary to fill template
    global_pf_keys = {
        'temp_step': temp_step,
        'ntemps': ntemps,
        'rel_temp_inc': rel_temp_inc,
        'atom_dist_min': atom_dist_min
    }

    # Set template name and path for the global keywords section for messpf run
    template_file_name = 'global_pf.mako'
    template_file_path = os.path.join(SECTION_PATH, template_file_name)

    # Build global section string
    global_pf_str = Template(filename=template_file_path).render(**global_pf_keys)

    return global_pf_str

"""
Writes the energy transfer section of a MESS input file
"""

import os
from mako.template import Template


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(SRC_PATH, 'templates')
SECTION_PATH = os.path.join(TEMPLATE_PATH, 'sections')

def write_energy_transfer(exp_factor, exp_power, exp_cutoff,
                          eps1, eps2,
                          sig1, sig2,
                          mass1, mass2):
    """ Writes the energy transfer section of the MESS input file
        :param float exp_factor: Exponent factor  
        :param float exp_power: Exponent power
        :param float exp_cutoff: Exponent cutoff  
        :param float eps1: Epsilon of Species 1 
        :param float eps2: Epsilon of Species 2
        :param float sig1: Sigma of Species 1 
        :param float sig2: Sigma of Species 2
        :param float mass1: Mass of Species 1 
        :param float mass2: Mass of Species 2
        :return energy_transfer_str: String for section
        :rtype: string
    """

    # Create dictionary to fill template
    energy_trans_keys = {
        'exp_factor': exp_factor,
        'exp_power': exp_power,
        'exp_cutoff': exp_cutoff,
        'epsilon1': eps1,
        'epsilon2': eps2,
        'sigma1': sig1,
        'sigma2': sig2,
        'mass1': mass1,
        'mass2': mass2
    }

    # Set template name and path for the energy transfer section 
    template_file_name = 'energy_transfer.mako'
    template_file_path = os.path.join(SECTION_PATH, template_file_name)

    # Build energy transfer string
    energy_trans_str = Template(filename=template_file_path).render(**energy_trans_keys)

    return energy_trans_str

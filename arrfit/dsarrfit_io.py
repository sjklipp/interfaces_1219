""" interface to arrfit code by sjk
"""

import os
import subprocess
from mako.template import Template


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))

def write_arrfit_inp(temps, ks, local_t):
    """ write the arrfit input file
    """

    # Format the lines with temps and rate constants
    tk_str = ''
    for i in range(len(temps)):
        tk_str += '{}  {}'.format(temps[i], ks[i])

    # Build the fill value dictionary
    fit_keys = {
        'local_T': local_t,
        'tk_lines': tk_str
        }

    # Set template name and path for the global keywords section
    template_file_name = 'dsarrfit.mako'
    template_file_path = os.path.join(SRC_PATH, template_file_name)

    # Build global section string
    arrfit_inp_str = Template(filename=template_file_path).render(**fit_keys)

    return arrfit_inp_str


def run_arrfit(path):
    """ run arrfit code
    """

    os.chdir(path)
    subprocess.check_call(['./dsarrfit.x_cfg'])


def parse_arrfit(output_string, conv_factor=1.987):
    """ obtain information from the arrfit output
    """

    # Get the lines in reverse order
    lines = output_string.readlines()
    lines.reverse()

    # Loop over the lines and find the resulting fit params
    for line in output_string.splitlines().reverse():
        if line.startswith(' params'):
            params1 = lines[lines.index(line)-1].split()
            params2 = lines[lines.index(line)-2].split()
            break

    # Combine the params and convert to floats
    best_guess = (
            [float(param) for param in params1] +
            [float(param) for param in params2]
    )

    # Convert units for final params based as needed
    best_guess[2] *= conv_factor
    best_guess[5] *= conv_factor

    return best_guess

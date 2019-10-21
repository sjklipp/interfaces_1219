""" interface to dsarrfit code by sjk
"""

import os
import re
import subprocess
from mako.template import Template


RC = 1.98720425864083e-3  # Gas Constant in kcal/mol.K

# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))


def write_input(temps, rate_constants,
                a_guess=8.1e-11,
                n_guess=-0.01,
                ea_guess=1000.0):
    """ write the dsarrfit input file
    """

    # Format the lines with temps and rate constants
    assert len(temps) == len(rate_constants)
    num_tk = len(temps)
    tk_str = ''
    for i, _ in enumerate(temps):
        tk_str += '{0:<10.1f}{1:<1.3E}\n'.format(temps[i], rate_constants[i])

    # Build the fill value dictionary
    fit_keys = {
        'a_guess': a_guess,
        'n_guess': n_guess,
        'ea_guess': ea_guess,
        'num_tk': num_tk,
        'tk_lines': tk_str
        }

    # Set template name and path for the global keywords section
    template_file_name = 'dsarrfit.mako'
    template_file_path = os.path.join(SRC_PATH, template_file_name)

    # Build global section string
    dsarrfit_str = Template(filename=template_file_path).render(**fit_keys)

    return dsarrfit_str


def run_dsarrfit(path):
    """ run arrfit code
    """
    start_path = os.getcwd()
    os.chdir(path)
    # Set the full path to the dsarrfit executable
    exe_path = os.path.join(SRC_PATH, 'dsarrfit', 'dsarrfit.x_cfg')
    # Run the executable
    subprocess.check_call([exe_path])
    os.chdir(start_path)


def read_params(output_string, fit, conv_factor=1.000):
    """ obtain information from the arrfit output
    """

    assert fit in ('single', 'double')

    # Loop over the lines and find the resulting fit params line
    lines = output_string.splitlines()
    lines.reverse()
    for line in lines:
        if line.startswith(' results for iteration'):
            params_str = lines[lines.index(line)-3]
            break

    # Grab the fitting parameters; multiply Ea param by given conversion factor
    # Multiple A by given conversion factor
    # Multiply Ea/R term by R to get Ea
    fit_params = [float(param) for param in params_str.split()]
    # print(fit_params)
    fit_params[0] *= conv_factor
    fit_params[2] *= RC
    if fit == 'double':
        fit_params[3] *= conv_factor
        fit_params[5] *= RC

    # print(fit_params)

    return fit_params


if __name__ == '__main__':
    with open('0/darrfit.out', 'r') as f:
        output_string = f.read()
    fit = 'double'
    read_params(output_string, fit, conv_factor=1.000)

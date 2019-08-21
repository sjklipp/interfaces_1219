""" interface to arrfit code by sjk
"""

import os
import subprocess
from mako.template import Template


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))

def write_arrfit_inp(temps, ks, locat_T):
    """ write the arrfit input file
    """

    # Format the lines with temps and rate constants
    tk_str = ''
    for i in range(len(temps)):
        tk_str += '{}  {}'.format(temps[i], ks[i])

    # Build the fill value dictionary
    arrfit_keys = {
        'local_T': locat_T,
        'tk_lines': tk_str
        }

    # Set template name and path for the global keywords section
    template_file_name = 'arrfit.mako'
    template_file_path = os.path.join(SRC_PATH, template_file_name)

    # Build global section string
    arrfit_inp_str = Template(filename=template_file_path).render(**arrfit_keys)

    return arrfit_inp_str


def run_arrfit(path):
    """ run arrfit code
    """

    os.chdir(path)
    subprocess.check_call(['./dsarrfit.x_cfg'])


def parse_arrfit(output_string):
    """ obtain information from the arrfit output
    """

    lines = output_string.readlines()
    lines.reverse()
    for line in lines:
        if line.startswith(' params'):
            bits1 = lines[lines.index(line)-1].split()
            bits2 = lines[lines.index(line)-2].split()
            break
    if (len(bits1) == 6) and (len(bits2) == 6):
        job = True
        best_guess = np.zeros(len(bits1), dtype=np.float64)
        for (b, bit) in enumerate(bits1):
            if '+' in bit:
                if bit[bit.index('+') - 1] != 'E':
                    bit = bit.replace('+', 'E+')
                # zador
                if len(bit) >= 9:
                    if (bit[8] == '-') and (bit[7] != 'E'):
                        bit = bit.replace('-', 'X', 1)
                        bit = bit.replace('-', 'E-')
                        bit = bit.replace('X', '-')
                if len(bit) >= 8:
                    if (bit[7] == '-') and (bit[6] != 'E'):
                        bit = bit.replace(bit[7], 'E-')

            best_guess[b] = float(bit)

        best_guess[2] = best_guess[2]*1.987
        best_guess[5] = best_guess[5]*1.987
    else:
        print('FAILED! Falling back to python solver')
        job = False


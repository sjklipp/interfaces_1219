#!/usr/bin/env python
"""
 Calculates heats-of-formation
"""

import os
import subprocess
from mako.template import Template

# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
TEMPLATE_PATH = os.path.dirname(os.path.realpath(__file__))


def write_thermp_input(formula, deltaH,
                       enthalpyT=0.0, breakT=1000.0,
                       thermp_file_name='thermp.dat'):
    """ Writes the input file for thermp
    """

    # Get the stoichiometry of all elements to build composition string
    composition_str = ''

    # Create a fill value dictionary
    thermp_keys = {
        'formula': formula,
        'deltaH': deltaH,
        'enthalpyT': enthalpyT,
        'breakT': breakT,
        'composition_str': composition_str
    }

    # Set template name and path for an atom
    template_file_name = 'thermp.mako'
    template_file_path = os.path.join(TEMPLATE_PATH, template_file_name)

    # Build a ProjRot input string
    thermp_str = Template(filename=template_file_path).render(**thermp_keys)

    # Write the file
    with open(thermp_file_name, 'w') as thermp_file:
        thermp_file.write(thermp_str)


def read_thermp_dHf(output_string):
    """
    Obtains deltaHf from thermp output
    """

    dHf298_str = ' h298 final\s*([\d,\-,\.]*)'
    dHf298 = float(re.findall(dHf298_str, output_string)[-1])

    return dHf298


def run_thermp(path, thermp_file_name='thermp.dat', pf_file_name='pf.out'):
    """
    Runs thermp.exe
    Requires pffile and thermpfile to be present
    """

    # Set full paths to files
    thermp_file = os.path.join(path, thermp_file_name)
    pf_file = os.path.join(path, pf_file_name)

    # Check for the existance of ThermP input and PF output
    assert os.path.exists(thermp_file)
    assert os.path.exists(pf_file)

    # Run thermp
    subprocess.check_call(['thermp.exe', thermp_file])


def run_pac99(path, formula):
    """
    Run pac99 for a given species name (formula)
    https://www.grc.nasa.gov/WWW/CEAWeb/readme_pac99.htm
    requires formula+'i97' and new.groups files
    """

    # Set file names for pac99
    i97_file = os.path.join(path, formula + '.i97')
    newgroups_file = os.path.join(path, 'newgroups')

    # Check for the existance of pac99 files
    assert os.path.exists(i97_file)
    assert os.path.exists(newgroups_file)

    # Run pac99
    subprocess.call(['pac99'])
    subprocess.call([formula])


def read_polyn(string, poly):
    """ read in the polyn
    """
    
    pac_polynomial = ''

    return pacc_polynomial


def convert_pac_to_chemkin(pac_polynomial_str):
    """ convert the polynimal from pac format to chemkin polynomial
    """

    # Get a reordered string
    polynomial_lines = pac_polynomial_str.splitlines()
    reorder_lines = [] 
    reorder_lines.append(polynomial_lines[9])
    reorder_lines.append(polynomial_lines[10])
    reorder_lines.append(polynomial_lines[6])
    reorder_lines.append(polynomial_lines[7])
    reorder_lines.append(polynomial_lines[3])
    reorder_lines.append(polynomial_lines[4])
    
    poly_vals = [val for line in reorder_lines
                     for val in line.split() 
                     if val != 0]

    chemkin_polynomial_str = ''
    nline = 2
    for i, val in enumerate(poly_vals):
        if i % 5 == 0:
            chemkin_polynomial_str += str(nline) + '\n'
            nline += 1
        else:
            chemkin_polynomial_str += str(val)

    return chemkin_polynomial

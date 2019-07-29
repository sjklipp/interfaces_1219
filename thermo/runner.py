"""
 Generates NASA Polynomial from MESS+THERMP+PAC99 outputs
"""

import os
import subprocess
from mako.template import Template
from . import util


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
TEMPLATE_PATH = os.path.dirname(os.path.realpath(__file__))


def write_thermp_input(formula, deltaH,
                       enthalpyT=0.0, breakT=1000.0,
                       thermp_file_name='thermp.dat'):
    """ Writes the input file for thermp
    """

    # Get the stoichiometry of all elements to build composition string
    atom_dict = util.get_atom_counts_dict(formula)
    composition_str = ''
    for key, val in atom_dict.items():
        composition_str += '{0}  {1}\n'.format(key, val)
    composition_str = composition_str.rstrip()

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


def run_thermp(path, thermp_file_name='thermp.dat', pf_file_name='pf.out'):
    """
    Runs thermp.exe
    Requires thermp input file to be present
    partition function (pf) output file and 
    """

    # Set full paths to files
    thermp_file = os.path.join(path, thermp_file_name)
    pf_file = os.path.join(path, pf_file_name)

    # Check for the existance of ThermP input and PF output
    assert os.path.exists(thermp_file)
    assert os.path.exists(pf_file)

    # Run thermp
    subprocess.check_call(['thermp', thermp_file])


def run_pac99(path, formula):
    """
    Run pac99 for a given species name (formula)
    https://www.grc.nasa.gov/WWW/CEAWeb/readme_pac99.htm
    requires formula+'i97' and new.groups files
    """

    # Set file names for pac99
    i97_file = os.path.join(path, formula + '.i97')
    newgroups_file = os.path.join(path, 'new.groups')

    # Check for the existance of pac99 files
    assert os.path.exists(i97_file)
    assert os.path.exists(newgroups_file)

    # Run pac99
    p = subprocess.Popen('pac99', stdin=subprocess.PIPE)
    p.communicate(formula)

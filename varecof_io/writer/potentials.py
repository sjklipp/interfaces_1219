"""
Writes the fortran files needed for the correction potential
"""

import os
import subprocess
from mako.template import Template
from mess_io.writer import util


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(SRC_PATH, 'templates')


def species(potentials, aidx, bidx, asym, bsym, species_name):
    """ Writes string for correction potential for some species Fortran file
        :return : String for the mol_corr.f file
        :rtype: string
    """

    assert potentials
    assert len(potentials[i]) same

    pot_label_str = ''
    if pot_labels:
        for i, label in enumerate(pot_labels):
            pot_label_str += 'c     dv{0} = {1} correction'.format(
                str(i+1), label)
    
    npot = len(potentials)
    npot_terms = len(potentials[0])
    dv_defs = ''
    for i in range(npot):
        dv_defs +='dv{0}({1})'.format(str(i+1), npot_terms) 
    
    rvals = ''
    for i, rval in enumerate(rvals):
        rvals += '     data rinp({0}) / {1} /'.format(
            str(i+1), rval)
    rmin = min(rvals)
    rmax = max(rvals)

    dv_vals = ''
    for i, potential in enumerate(potentials):
        for j, term in enumerate(potential):
            dv_vals += '     data dv{0}({1}) / {2} /'.format(
                str(i+1), str(j+1), term)

    spline_str = ''
    for i in range(npot): 
        if i == 0:
            splint_str += '      if (ipot.eq.{0}) then'.format(str(i+1))) 
        else:
            splint_str += '      else if (ipot.eq.{0}) then'.format(str(i+1))) 
        spline_str += '      call spline(rinp,dv{0},nrin,dvp1,dvpn,dv20)'.format(
            str(i+1))
        spline_str += '      call splint(rinp,dv{0},dv20,nrin,r{1}{2},{3})'.format(
            str(i+1), asym, bsym, species_name+'_corr'))
    spline_str += 'endif'

    # Create dictionary to fill template
    corr_keys = {
        'species_name': species_name,
        'asym': asym,
        'bsym': bsym,
        'corr_labels_str': corr_labels_str,
        'npot': npot,
        'npot_terms': npot_terms,
        'dv_defs': dv_defs,
        'rvals': rvals,
        'dv_vals': dv_vals,
        'rmin': rmin,
        'rmax': rmax,
        'aidx': aidx,
        'bidx': bidx,
        'spline_str': spline_str
    }

    # Set template name and path for the mol_corr template
    template_file_name = 'mol_corr.mako'
    template_file_path = os.path.join(SPECIES_PATH, template_file_name)

    # Build mol_corr.f string
    mol_corr_str = Template(filename=template_file_path).render(**corr_keys)

    return mol_corr_str


def dummy():
    """ Writes string for the dummy correction potention potential Fortran file
        :return : String for the dummy_corr.f file
        :rtype: string
    """

    # Set template name and path for the dummy_corr template
    template_file_name = 'dummy_corr.mako'
    template_file_path = os.path.join(SPECIES_PATH, template_file_name)

    # Build dummy_corr.f string
    dummy_corr_str = Template(filename=template_file_path).render()

    return dummy_corr_str


def auxiliary():
    """ Writes string for the potential auxiliary functions Fortran file
        :return : String for the pot_aux.f file
        :rtype: string
    """

    # Set template name and path for the pot_aux template
    template_file_name = 'pot_aux.mako'
    template_file_path = os.path.join(SPECIES_PATH, template_file_name)

    # Build pot_aux.f string
    pot_aux_str = Template(filename=template_file_path).render()

    return pot_aux_str


def makefile():
    """ Writes string for a makefile to compile correction potentials
        :return : String for the makefile
        :rtype: string
    """

    # Set template name and path for the pot_aux template
    template_file_name = 'makefile.mako'
    template_file_path = os.path.join(SPECIES_PATH, template_file_name)

    # Build pot_aux.f string
    makefile_str = Template(filename=template_file_path).render()

    return makefile_str


def compile_correction_pot(path):
    """ compile the correction potential
    """
    subprocess(['make'], cwd=path)

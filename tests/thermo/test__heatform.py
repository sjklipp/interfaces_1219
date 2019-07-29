"""
Tests calculating the 0 K heat-of-formation
"""

import os
from thermo import heatform
from thermo import util


# Inchi string for methyl nitrate (CH3ONO2)
ICH = 'InCHI=1S/CH3NO3/c1-5-2(3)4/h1H3'

# Thermp output file name
THERMP_OUTFILE_NAME = os.path.join(os.getcwd(), 'run', 'thermp.out')


def test__calc_hform_0k():
    """ calculates 0 K heat-of-formation for a species
    """

    # Get the molecular formula from the inchi string
    formula = util.inchi_formula(ICH)
    print('\nformula:')
    print(formula)

    # Get atom count dictionary
    atom_dict = util.get_atom_counts_dict(formula)
    print('\natom dict:')
    print(atom_dict)

    # Get the list of the basis
    basis = heatform.select_basis(atom_dict)
    print('\nbasis:')
    print(basis)

    # Get the coefficients for the balanced heat-of-formation eqn
    coeff = heatform.calc_coefficients(basis, atom_dict)
    print('\ncoeff:')
    print(coeff)

    # Get the energy for the species
    e_mol = heatform.get_mol_energy()
    print('\ne_mol:')
    print(e_mol)

    # Get the energy for each of the basis species
    e_basis = heatform.get_basis_energy(basis)
    print('\ne_basis:')
    print(e_basis)

    # Get the 0 K heat of formation
    hform = heatform.calc_hform_0k(e_mol, e_basis, coeff)
    print('\nhform(0 K):')
    print(hform)


def test__read_hform_298k():
    """ reads the 298 K heat-of-formation value from thermp output
    """

    # Read the thermp output
    with open(THERMP_OUTFILE_NAME, 'r') as thermp_outfile:
        thermp_out_str = thermp_outfile.read()

    # Get the 0 K heat of formation
    hform = heatform.get_hform_298k_thermp(thermp_out_str)
    print('\nhform(298 K):')
    print(hform)


if __name__ == '__main__':
    test__calc_hform_0k()
    test__read_hform_298k()

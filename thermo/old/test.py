"""
Driver routine for calculating the heat-of-formation at 0 K
"""

import util
import heatform

# Set init variables
ICH = 'InCHI=1S/CH3NO3/c1-5-2(3)4/h1H3'
FORMULA = util.inchi_formula(ICH)
print('\nformula:')
print(FORMULA)

# Get atom count dictionary
ATOM_DCT = util.get_atom_counts_dict(FORMULA)
print('\natom dict:')
print(ATOM_DCT)

# Get the list of the basis
BASIS = heatform.select_basis(ATOM_DCT)
print('\nbasis:')
print(BASIS)

# Get the coefficients for the balanced heat-of-formation eqn
COEFF = heatform.calc_coefficients(BASIS, ATOM_DCT)
print('\ncoeff:')
print(COEFF)

# Get the energy for the species
E_MOL = heatform.get_mol_energy()
print('\ne_mol:')
print(E_MOL)

# Get the energy for each of the basis species
E_BASIS = heatform.get_basis_energy(BASIS)
print('\ne_basis:')
print(E_BASIS)

# Get the 0 K heat of formation
HFORM = heatform.calc_hform_0k(E_MOL, E_BASIS, COEFF)
print('\nhform:')
print(HFORM)

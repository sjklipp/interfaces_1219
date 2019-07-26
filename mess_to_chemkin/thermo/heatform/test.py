"""
Driver routine for thermo
"""

import util
import heatform

# import automol.inchi
# ICH = 'InChI=1/C2H6O/c1-2-3/h3H,2H2,1H3'
# FORMULA = automol.inchi.formula_sublayer(ICH)
# print(FORMULA)

# Set init variables
FORMULA = 'C2H6O'

# Get atom count dictionary
ATOM_DCT = util.get_atom_counts_dict(FORMULA)
print(ATOM_DCT)

# Get the list of the basis
BASIS = heatform.select_basis(ATOM_DCT)
print(BASIS)

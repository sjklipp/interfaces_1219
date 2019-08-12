"""
Computes the Heat of Formation at 0 K for a given species
"""

import os
import csv
import numpy as np
from qcelemental import constants as qcc
import autoparse.pattern as app
import autoparse.find as apf
from . import util


KJ2KCAL = qcc.conversion_factor('kJ/mol', 'kcal/mol')

SRC_PATH = os.path.dirname(os.path.realpath(__file__))


def get_hform_298k_thermp(output_string):
    """
    Obtains deltaHf from thermp output
    """

    dHf298_pattern = ('h298 final' +
                  app.one_or_more(app.SPACE) +
                  app.capturing(app.FLOAT))
    dHf298 = float(apf.last_capture(dHf298_pattern, output_string))

    return dHf298


def calc_hform_0k(hzero_mol, hzero_basis, basis, coeff, ref_set):
    """ calculates the heat-of-formation at 0 K
    """

    # Calculate the heat of formation
    dhzero = hzero_mol
    for i, spc in enumerate(basis):
        h_basis = get_ref_h(spc, ref_set, 0)
        if h_basis is None:
            h_basis = 0.0
        dhzero += coeff[i] * h_basis * KJ2KCAL
        dhzero -= coeff[i] * hzero_basis[i]

    return dhzero


def get_ref_h(species, ref, temp):
    """ gets a reference value
    """

    # Set path and name to thermo database file
    thermodb_name = 'thermodb_{}K.csv'.format(str(int(temp)))
    thermodb_file = os.path.join(SRC_PATH, thermodb_name)

    # Find the energy value for the given species and enery type
    h_species = None
    with open(thermodb_file, 'r') as db_file:
        reader = csv.DictReader(db_file)
        for row in reader:
            if row['inchi'] == species:
                val = row[ref]
                if val == '':
                    val = None
                h_species = float(val)

    return h_species


def select_basis(atom_dct, attempt=0):
    """
    Given a list of atoms, generates a list of molecules
    that is best suited to serve as a basis for those atoms

    :param atomlist: list of atoms
    :type atomlist: list
    :param attempt: ???
    :type attempt: ???

    OUPUT:
    basis    - recommended basis as a list of stoichiometries
    """

    # Determine number of basis species required
    nbasis = len(atom_dct)

    # Get a list of all the atom types in the molecule
    atoms = list(atom_dct.keys())

    # Create list of inchi keys corresponding to basis species
    basis = []
    counter = 1
    # N2
    if 'N' in atoms and attempt < 2 and counter <= nbasis:
        basis.append('InChI=1S/N2/c1-2')
        counter += 1
    # NH3
    if 'N' in atoms and 'H' in atoms and attempt > 1 and counter <= nbasis:
        basis.append('InChI=1S/H3N/h1H3')
        counter += 1
    # SO2
    if 'S' in atoms and counter <= nbasis:
        basis.append('InChI=1S/O2S/c1-3-2')
        counter += 1
    # H2
    if 'H' in atoms and attempt < 2 and counter <= nbasis:
        basis.append('InChI=1S/H2/h1H')
        counter += 1
    # H2
    elif 'H' in atoms and 'C' not in atoms and attempt < 3 and counter <= nbasis:
        basis.append('InChI=1S/H2/h1H')
        counter += 1
    # O2
    if 'O' in atoms and attempt < 3 and counter <= nbasis:
        basis.append('InChI=1S/O2/c1-2')
        counter += 1
    # CH4
    if 'C' in atoms and attempt < 4 and counter <= nbasis:
        basis.append('InChI=1S/CH4/h1H4')
        counter += 1
    # H2O
    if 'O' in atoms and 'H' in atoms and attempt < 4 and counter <= nbasis:
        basis.append('InChI=1S/H2O/h1H2')
        counter += 1
    # CO2
    if 'C' in atoms and 'O' in atoms and attempt < 5 and counter <= nbasis:
        basis.append('InChI=1S/CO2/c2-1-3')
        counter += 1
    # CH2O
    if 'C' in atoms and 'O' in atoms and attempt < 5 and counter <= nbasis:
        basis.append('InChI=1S/CH2O/c1-2/h1H2')
        counter += 1
    # CH3OH
    if 'C' in atoms and 'O' in atoms and counter <= nbasis:
        basis.append('InChI=1S/CH4O/c1-2/h2H,1H3')
        counter += 1
    # CH3CH3
    if 'C' in atoms and counter <= nbasis:
        basis.append('InChI=1S/C2H6/c1-2/h1-2H3')
        counter += 1
    # SO2
    if 'S' in atoms and counter <= nbasis:
        basis.append('InChI=1S/O2S/c1-3-2')
        counter += 1
    # H2
    if 'H' in atoms and attempt < 1 and counter <= nbasis:
        basis.append('InChI=1S/H2/h1H')
        counter += 1
    # H2
    elif 'H' in atoms and 'C' not in atoms and attempt < 3 and counter <= nbasis:
        basis.append('InChI=1S/H2/h1H')
        counter += 1
    # O2
    if 'O' in atoms and attempt < 2 and counter <= nbasis:
        basis.append('InChI=1S/O2/c1-2')
        counter += 1
    # CH4
    if 'C' in atoms and attempt < 3 and counter <= nbasis:
        basis.append('InChI=1S/CH4/h1H4')
        counter += 1
    # H2O
    if 'O' in atoms and 'H' in atoms and attempt < 3 and counter <= nbasis:
        basis.append('InChI=1S/H2O/h1H2')
        counter += 1
    # CO2
    if 'C' in atoms and 'O' in atoms and attempt < 4 and counter <= nbasis:
        basis.append('InChI=1S/CO2/c2-1-3')
        counter += 1
    # CH2O
    if 'C' in atoms and 'O' in atoms and attempt < 4 and counter <= nbasis:
        basis.append('InChI=1S/CH2O/c1-2/h1H2')
        counter += 1
    # CH3OH
    if 'C' in atoms and 'O' in atoms and counter <= nbasis:
        basis.append('InChI=1S/CH4O/c1-2/h2H,1H3')
        counter += 1

    return basis


def get_reduced_basis(basis_formulae, species_formula):
    """
    Form a matrix for a given basis and atomlist
    INPUT:
    input_basis     - ich stringes for set of reference molecules
    atomlist  - list of atoms (all atoms that appear
                in basis should be in atomlist)
    OUTPUT:
    mat       - matrix (length of basis by length of atomlist)
                (square if done right)
    """
    
    # Get the basis formulae list
    #basis_formulae = [util.inchi_formula(spc) for spc in basis]

    reduced_basis = []
    for i, basis_formula in enumerate(basis_formulae):
        basis_atom_dict = util.get_atom_counts_dict(basis_formula)
        flag = True
        for key, _ in basis_atom_dict.items():
            if key not in species_formula:
                flag = False

        if flag:
            reduced_basis.append(basis_formulae[i])

    return reduced_basis


def calc_coefficients(basis, mol_atom_dict):
    """
    Form a matrix for a given basis and atomlist
    INPUT:
    basis     - basis of molecules
    atomlist  - list of atoms (all atoms that appear
                in basis should be in atomlist)
    OUTPUT:
    mat       - matrix (length of basis by length of atomlist)
                (square if done right)
    """

    # Initialize an natoms x natoms matrix
    nbasis = len(basis)
    basis_mat = np.zeros((nbasis, nbasis))

    # Get the basis formulae list
    basis_formulae = [util.inchi_formula(spc) for spc in basis]

    # Set the elements of the matrix
    for i, spc in enumerate(basis_formulae):
        basis_atom_dict = util.get_atom_counts_dict(spc)
        basis_vals = []
        for key in mol_atom_dict.keys():
            if key in basis_atom_dict:
                basis_vals.append(basis_atom_dict[key])
            else:
                basis_vals.append(0)
        basis_mat[i] = basis_vals

    #  Transpose
    basis_mat = basis_mat.T

    # Form stoich vector
    stoich_vec = np.zeros(len(mol_atom_dict))
    for i, key in enumerate(mol_atom_dict.keys()):
        stoich_vec[i] = mol_atom_dict[key]

    # Solve C = M^-1 S
    basis_mat = np.linalg.inv(basis_mat)
    coeff = np.dot(basis_mat, stoich_vec)

    return coeff

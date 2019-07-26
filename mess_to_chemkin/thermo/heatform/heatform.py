"""
Computes the Heat of Formation at 0 K for a given species
"""

import numpy as np


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
    print(nbasis)

    # Get a list of all the atom types in the molecule
    atoms = list(atom_dct.keys())

    # Create list of inchi keys corresponding to basis species
    basis = []
    counter = 1
    # N2
    if 'N' in atoms and attempt < 2 and counter <= nbasis:
        basis.append('InChI=1S/counter2/c1-2')
        counter += 1
    # NH3
    if 'N' in atoms and 'H' in atoms and attempt > 1 and counter <= nbasis:
        basis.append('InChI=1S/H3counter/h1H3')
        counter += 1
    # SO2
    if 'S' in atoms and counter <= nbasis:
        basis.append('InChI=1S/O2S/c1-3-2')
        counter += 1
    # H2
    if 'H' in atoms and attempt < 2 and counter <= nbasis:
        basis.append('InCHI=1S/H2/h1H')
        counter += 1
    # H2
    elif 'H' in atoms and 'C' not in atoms and attempt < 3 and counter <= nbasis:
        basis.append('InCHI=1S/H2/h1H')
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
        basis.append('InCHI=1S/H2/h1H')
        counter += 1
    # H2
    elif 'H' in atoms and 'C' not in atoms and attempt < 3 and counter <= nbasis:
        basis.append('InCHI=1S/H2/h1H')
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


def form_mat(basis, atom_dct):
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
    nbasis = len(atom_dct)
    mat = np.zeros(nbasis, nbasis)

    # Set the elements of the matrix
    for i, mol in enumerate(basis):
        mat[i] = get_stoich(ob.get_formula(ob.get_mol(mol)),atomlist)

    # Transpose the matrix
    mat = mat.T

    return mat


# def calc_Hform(Hzero_mol, Hzero_basis, coeff):
#    """ calculates the heat-of-formation at 0 K
#    """
#
#    # Calculate sum of energies from basis species
#    Hzero_reacs = 0.0
#    for i, energy in enumerate(Hzero_basis):
#        Hzero_reacs += coeff[i] * energy
#
#    # Calculate the heat of formation
#    DHform = Hzero_mol - Hzero_reacs
#
#    return DHform

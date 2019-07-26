import re
import numpy as np
import os
import sys
from . import iotools as io
from . import qctools as qc
from . import patools as pa
from . import obtools as ob
from . import unittools as ut
import logging


def get_atomlist(mol):
    """
    Makes a list of all atoms in a molecule
    INPUT:
    mol      - stoichiometry of molecule
    OUTPUT:
    atomlist - list of distinct atoms in that molecule
    """
    atomlist = []
    elements = {'He','Li','Be','Ne','Na','Mg','Al','Si','Cl','Ar',
      'Ca','Sc','Ti','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge',
      'As','Se','Br','Kr','C','B','H','O','F','S','N','P','K','V'}
    for el in elements:
        if el in mol:
           atomlist.append(el)
           mol = mol.replace(el,'')
    return atomlist


def select_basis(atomlist, attempt=0):
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
    count = len(atomlist)-1
    basis = []
    i = 0
    # N2
    if 'N' in atomlist and i <= count and attempt < 2:
        basis.append('InChI=1S/N2/c1-2')
        i += 1
    # NH3
    if 'N' in atomlist and 'H' in atomlist and i <= count and attempt > 1:
        basis.append('InChI=1S/H3N/h1H3')
        i += 1
    # SO2
    if 'S' in atomlist and i <= count:
        basis.append('InChI=1S/O2S/c1-3-2')
        i += 1
    # H2
    if 'H' in atomlist and i <= count and attempt < 2:
        basis.append('InCHI=1S/H2/h1H')
        i += 1
    # H2
    elif 'H' in atomlist and 'C' not in atomlist and i <= count and attempt < 3:
        basis.append('InCHI=1S/H2/h1H')
        i += 1
    # O2
    if 'O' in atomlist and i <= count and attempt < 3:
        basis.append('InChI=1S/O2/c1-2')
        i += 1
    # CH4
    if 'C' in atomlist and i <= count and attempt < 4:
        basis.append('InChI=1S/CH4/h1H4')
        i += 1
    # H2O
    if 'O' in atomlist and 'H' in atomlist and i <= count and attempt < 4:
        basis.append('InChI=1S/H2O/h1H2')
        i += 1
    # CO2
    if 'C' in atomlist and 'O' in atomlist and i <= count and attempt < 5:
        basis.append('InChI=1S/CO2/c2-1-3')
        i += 1
    # CH2O
    if 'C' in atomlist and 'O' in atomlist and i <= count and attempt < 5:
        basis.append('InChI=1S/CH2O/c1-2/h1H2')
        i += 1
    # CH3OH
    if 'C' in atomlist and 'O' in atomlist and i <= count:
        basis.append('InChI=1S/CH4O/c1-2/h2H,1H3')
        i += 1
    # CH3CH3
    if 'C' in atomlist and i <= count:
        basis.append('InChI=1S/C2H6/c1-2/h1-2H3')
        i += 1
    # SO2
    if 'S' in atomlist and i <= count:
        basis.append('InChI=1S/O2S/c1-3-2')
        i += 1
    # H2
    if 'H' in atomlist and i <= count and attempt < 1:
        basis.append('InCHI=1S/H2/h1H')
        i += 1
    # H2
    elif 'H' in atomlist and 'C' not in atomlist and i <= count and attempt < 3:
        basis.append('InCHI=1S/H2/h1H')
        i += 1
    # O2
    if 'O' in atomlist and i <= count and attempt < 2:
        basis.append('InChI=1S/O2/c1-2')
        i += 1
    # CH4
    if 'C' in atomlist and i <= count and attempt < 3:
        basis.append('InChI=1S/CH4/h1H4')
        i += 1
    # H2O
    if 'O' in atomlist and 'H' in atomlist and i <= count and attempt < 3:
        basis.append('InChI=1S/H2O/h1H2')
        i += 1
    # CO2
    if 'C' in atomlist and 'O' in atomlist and i <= count and attempt < 4:
        basis.append('InChI=1S/CO2/c2-1-3')
        i += 1
    # CH2O
    if 'C' in atomlist and 'O' in atomlist and i <= count and attempt < 4:
        basis.append('InChI=1S/CH2O/c1-2/h1H2')
        i += 1
    # CH3OH
    if 'C' in atomlist and 'O' in atomlist and i <= count:
        basis.append('InChI=1S/CH4O/c1-2/h2H,1H3')
        i += 1
    # CH3CH3
    if 'C' in atomlist and i <= count:
        basis.append('InChI=1S/C2H6/c1-2/h1-2H3')
        i += 1

    return basis


def get_stoich(mol, atomlist):
    """
    Given a molecule's stoichiometry and a list of atoms, finds the
    number of each atom that the molecule contains
    INPUT:
    mol       - molecule stoichiometry
    atomlist  - list of atoms

    OUTPUT:
    stoich    - list of numbers corresponding to the number of each atom
                in the atomlist that the molecule contains
    """
    stoichlist = np.zeros(len(atomlist))

    for i, atom in enumerate(atomlist):
        val = 0
        if atom in mol:
            a = re.compile(atom + '(\d*)')
            a = a.findall(mol)
            for b in a:
                if b == '':
                    b = '1'
                val += float(b)
        stoichlist[i] = val

    return stoichlist


def form_mat(basis, atomlist):
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
    mat = np.zeros((len(atomlist),len(atomlist)))
    for i,mol in enumerate(basis):
        mat[i] = get_stoich(ob.get_formula(ob.get_mol(mol)),atomlist)
    mat = mat.T

    return mat


def comp_coeff(mat, stoich):
    """
    Finds the coefficients that solve C = M^-1 S:
    C are the coefficients [a, b, c, d..] for
        basis [R1, R2, R3, R4...] that gives
        CxHyOzNw =  aR1 + bR2 + cR3 + dR4, where
        [x, y, z, w...] is S, our stoichiometry vector for a molecule.
    M is a nonsingular matrix that puts the basis in terms of a list of atoms
    param: mat: nonsingular matrix of ???
    type: mat: ???
    param: stoich: numbers giving stoichiometry in terms of a list of atoms
    type: stoich: list of ints
    OUTPUT:
    coeff  - coefficients [a,b,c,d...] as described above
    """

    mati = np.linalg.inv(mat)
    coeff = np.dot(mati, stoich)

    return coeff


def E_from_hfbasis(mol, basis, coefflist, E, opt, en, freq, anharm,dbdir='./'):
    """
    Uses the coefficients [a,b,c...] obtained from C = M^-1 S to find
    delH(CxHyOz) = adelH(R1) + bdelH(R2) + cdelH(R3) + Eo(CxHyOz) - aEo(R1) - bEo(R2) -cEo(R3)
    where Rn are our basis molecules, delH(Rn) are their heats of formation, and Eo(Rn) are their
    electronic energies computed at the same level of theory as Eo(CxHyOz)
    INPUTS:
    mol       - molecule named stoichiometrically
    basis     - selected basis molecule list
    coefflist - coefficients [a,b,c,d...] described above
    E         - electronic energy of molecule
    OUTPUTS:
    E        - 0K heat of formation of molecule

    """
    for i, bas in enumerate(basis):
        h   = nest_2_dic(bas,'delHf',  0)
        if h is None:
            h = 0
        E  +=  coefflist[i] * h * ut.kj2au
        e    =  find_E(bas, opt, en, freq, anharm=anharm,dbdir=dbdir)
        E   -=  coefflist[i] * e

    return E


def E_hfbasis_QTC(mol,basis,coefflist,E,opt, en, freq, parameters):
    """
    Uses the coefficients [a,b,c...] obtained from C = M^-1 S to find
    delH(CxHyOz) = adelH(R1) + bdelH(R2) + cdelH(R3) + Eo(CxHyOz) - aEo(R1) - bEo(R2) -cEo(R3)
    where Rn are our basis molecules, delH(Rn) are their heats of formation, and Eo(Rn) are their
    electronic energies computed at the same level of theory as Eo(CxHyOz)
    INPUTS:
    mol       - molecule named stoichiometrically
    basis     - selected basis molecule list
    coefflist - coefficients [a,b,c,d...] described above
    E         - electronic energy of molecule
    OUTPUTS:
    E        - 0K heat of formation of molecule

    """
    for i,bas in enumerate(basis):
        bas = ob.get_slabel(bas)
        smilesname = ob.get_smiles_filename(bas)
        formula   = ob.get_formula(bas)
        bas = bas.split('_')[0]
        h  =  nest_2_dic(bas,'delHf',  0)
        if h is None:
            smilesdir = io.join_path(parameters['database'], formula, smilesname)
            import os
            rundir =  io.join_path(*[os.getcwd().split(parameters['database'])[0],smilesdir,parameters['qcdirectory']])
            hoffile = io.join_path(*[rundir, formula + '.hofk'])
            if io.check_file(hoffile):
                h = float(io.read_file(hoffile).splitlines()[2].split()[0]) / ut.au2kcal / ut.kj2au
                logging.info('{} {}: {:5f}  pulled from {}'.format(bas, 'delHf', h, rundir))
            else:
                h = 0
        E  +=  coefflist[i] * h * ut.kj2au
        e    =  E_QTC(bas, opt, en, freq, parameters)
        if parameters['bac']:
            e += E_BAC(bas, parameters) / ut.au2kcal
        E   -=  coefflist[i] * e
    return E


def check(clist, basis,stoich,atomlist):
    """
    Makes sure nothing funky happened while computing coefficients
    """
    check = np.zeros(len(clist))
    statement = 'Coefficients produce correct stoichiometry\n'
    for i, c in enumerate(clist):
       check += c * get_stoich(ob.get_formula(ob.get_mol(basis[i])),atomlist)
    for i, sto in enumerate(stoich):
        if not check[i] == sto:
            statement = 'Coefficients do NOT produce correct stoichiometry'
            break
    return statement


def comp_coefficients(molform, basis='auto'):
    """ Calculates the coefficients for calculating the heat-of-formation
    """

    # Basis selection
    basprint = 'manually select basis'
    atomlist = get_atomlist(molform)
    basisselection = 0
    if is_auto(basis[0]):
        basis = select_basis(atomlist)
        basisselection += 1
        basprint = 'automatically generate basis'
    elif basis[0] == 'basis.dat':
        basis = io.read_file('basis.dat').split()
        basprint = 'read basis from basis.dat'

    for bas in basis:
        bas = ob.get_formula(ob.get_mol(bas))
        atomlist.extend(get_atomlist(   bas))

    # Compute atomlist, stoichlist, matrix, and coefficients
    atomlist = list(set(atomlist))
    stoich = get_stoich(molform,atomlist)
    mat = form_mat(basis,atomlist)

    # Pick a new basis if current one produces singular matrix
    for i in range(5):
        if np.linalg.det(mat) != 0:
            break

        basprint += '\nMatrix is singular -- select new basis'

        atomlist = get_atomlist(molform)
        if 'H' not in atomlist:
            atomlist.append('H')
        basis = select_basis(atomlist,basisselection)
        basisselection += 1

        for bas in basis:
            bas = ob.get_formula(ob.get_mol(bas))
            atomlist.extend(get_atomlist(bas))

        atomlist = list(set(atomlist))
        stoich   = get_stoich(molform,atomlist)
        mat      = form_mat(basis,atomlist)
        basprint +='\n\nBasis is: ' + ', '.join(basis)
        logging.debug( basprint)
        #basprint +='\n'.join(['\t'.join([{}.format(el) for el in line] for line in mat])
    basprint += '\n  ' + molform + '\t\t'
    basprint += '\t'.join([ob.get_formula(bas) for bas in basis])
    for i in range(len(mat)):
       basprint += '\n' + atomlist[i] + '  '
       basprint += str(stoich[i]) + '    \t'
       for el in mat[i]:
           basprint += str(el) + '\t'

    clist =  comp_coeff(mat,stoich)

    return clist, basis, basprint


def E_BAC(bas, parameters):
    ### Check dictionary ###
    from .heatform_db import db
    from . import qctools as qc
    from . import iotools as io
    slabel = qc.get_slabel(bas)
    calcindex = parameters['calcindex']
    qckeyword = parameters['qckeyword']
    qlabel = qc.get_qlabel(qckeyword, calcindex)
    bac = 0.
    if 'bac' in parameters['all results'][slabel][qlabel]:
        bac = parameters['all results'][slabel][qlabel]['bac']
    if bac:
        logging.debug('BAC for {0} {1} = {2} kcal/mol'.format(slabel, qlabel,bac))
    else:
        logging.error('BAC not found for {0} {1}'.format(slabel, qlabel))
    return  float(bac)


def E_QTC(bas, opt, en, freq, parameters):
    ### Check dictionary ###
    from .heatform_db import db
    from . import qctools as qc
    from . import iotools as io
    natom = ob.get_natom(bas)
    slabel = qc.get_slabel(bas)
    parameters['natom'] = natom
    calcindex = parameters['calcindex']
    qckeyword = parameters['qckeyword']
    qlabel = qc.get_qlabel(qckeyword, calcindex)
    en, zpve, bac = 0., 0., 0.

    if qlabel in parameters['all results'][slabel]:
        if 'energy' in parameters['all results'][slabel][qlabel]:
            en = parameters['all results'][slabel][qlabel]['energy']
    if en:
        logging.debug('Energy for {0} {1} = {2} Hartree'.format(slabel, qlabel,en))
    else:
        logging.error('Energy not found for {0} {1}'.format(slabel, qlabel))
    if qlabel in parameters['all results'][slabel]:
        if 'azpve' in parameters['all results'][slabel][qlabel]:
            zpve = parameters['all results'][slabel][qlabel]['azpve']
            zpvelabel = 'anharmonic ' + qlabel
        elif 'zpve' in parameters['all results'][slabel][qlabel]:
            zpve = parameters['all results'][slabel][qlabel]['zpve']
            zpvelabel = 'harmonic ' + qlabel
    else:
        for i in range(calcindex):
            qlabel = qc.get_qlabel(qckeyword, i)
            if qlabel in parameters['all results'][slabel]:
                if 'azpve' in parameters['all results'][slabel][qlabel]:
                    zpve = parameters['all results'][slabel][qlabel]['azpve']
                    zpvelabel = 'anharmonic ' + qlabel
                elif 'zpve' in parameters['all results'][slabel][qlabel]:
                    zpve = parameters['all results'][slabel][qlabel]['zpve']
                    zpvelabel = 'harmonic ' + qlabel
    if zpve:
        logging.debug('ZPVE (harmonic) for {0} {1} = {2} Hartree'.format(slabel,zpvelabel,zpve))
    else:
        logging.warning('ZPVE not found for {0} {1}'.format(slabel, qlabel))
    return  float(en) + float(zpve)

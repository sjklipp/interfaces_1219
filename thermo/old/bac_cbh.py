###CBH
def get_balance(smiles, frags):
    stoichs = {}
    for frag in frags:
       stoich = qc.get_atom_stoich(ob.get_formula(frag))
       for atom in stoich:
           if atom in stoichs:
               stoichs[atom] += stoich[atom] * frags[frag]
           else:
               stoichs[atom]  = stoich[atom] * frags[frag]
    stoich = qc.get_atom_stoich(ob.get_formula(smiles))
    balance = {}
    for atom in stoich:
        if atom in stoichs:
            balance[atom] = stoich[atom] - stoichs[atom]
        else:
            balance[atom] = stoich[atom] 
    return balance

def make_balanced(smiles, frags):
    balance = get_balance(smiles, frags)
    if 'C' in balance:
        if balance['C'] != 0:
            if not 'C' in frags:
                frags['C'] = balance['C'] 
            else:
                frags['C'] += balance['C'] 
    if 'N' in balance:
        if balance['N'] != 0:
            if not 'N' in frags:
                frags['N'] = balance['N'] 
            else:
                frags['N'] += balance['N'] 
    if 'O' in balance:
        if balance['O'] != 0:
            if not 'O' in frags:
                frags['O'] = balance['O'] 
            else:
                frags['O'] = balance['O']

    balance = get_balance(smiles, frags)
    if 'H' in balance:
        if balance['H'] != 0:
            if not '[H][H]' in frags:
                frags['[H][H]'] = balance['H'] / 2.
            else:
                frags['[H][H]'] += balance['H'] / 2.
    return frags

def isUnbalanced(balance):
    bal = True
    for atom in balance:
       if abs(balance[atom]) > 0:
           bal = False
    if bal:
       return False
    return True

def get_adj_mat(lines):
    lines = lines.lower().replace('\n\n','\n')
    lines = lines.splitlines()
    rows = []
    rad  = []
    rads = []
    startkey = 'molecular structure:'
    endkey = 'z-matrix atom order:'
    endkey2 = 'free radical'
    save = False
    for line in lines:
        if line.startswith(endkey):
            save = False
        if line.startswith(endkey2):
            save = False
            if 'is' in line:
                rad = line.split('site is ')[1]
            elif 'are' in line:
                rad = line.split('sites are ')[1].split()
        if save:
            rows.append(line)
        if line.startswith(startkey):
            save = True
    
    if rad:
        for i, atom in enumerate(rows[0].split()[1:]):
            if atom.strip() in rad:
                rads.append(i)
    atoms = qc.get_atom_list(rows[0])
    mat = []
    for row in rows[1:]:
        if row:
            mat.append(row.split()[1:-1])
    return atoms, mat, rads

def bond_sum(mat, atoms, heavy, i):
    bonds = 0
    for j, col in enumerate(mat[i]):
        if atoms[j] in heavy:
            if float(col) > 0:
                bonds += int(round(float(col)-0.5))
    return bonds

def lhs_rhs(frags):
    rhs = {}
    lhs = {}
    for frag in frags:
        if frags[frag] > 0:
            rhs[frag] = frags[frag]
        elif frags[frag] < 0:
            lhs[frag] = - frags[frag]
    return lhs, rhs

def print_lhs_rhs(smiles, frags):
    lhs, rhs = lhs_rhs(frags)
    lhsprint = ob.get_formula(smiles)
    rhsprint = ''
    for frag in rhs:
        if rhsprint:
            rhsprint += ' +  {:.1f} {} '.format( rhs[frag], ob.get_formula(frag))
        else:
            rhsprint  = ' {:.1f} {} '.format( rhs[frag], ob.get_formula(frag))
    for frag in lhs:
            lhsprint += ' +  {:.1f} {} '.format( lhs[frag], ob.get_formula(frag))
    return '{} --> {}'.format(lhsprint, rhsprint)

def CBHzed(smiles, mat, atoms, heavy, bond, rads):
    frags = {}
    valence = {'C':4, 'O':2, 'N':3}
    for i, row in enumerate(mat):
       if atoms[i] in heavy:
          val = valence[atoms[i]]
          frag = ''
          frag = atoms[i]
          if i in rads:
              val -= 1
          saturate = ''
          if val == 1:
              saturate = 'H'
          elif val > 1:
              saturate = 'H{:g}'.format(val)
          frag = '[{}{}]'.format(frag, saturate)
          frag = ob.get_slabel(frag).split('_')[0]
          if frag in frags:
              frags[frag] += 1
          else: 
              frags[frag] = 1
    
    frags = make_balanced(smiles, frags)
    return frags

def CBHone(smiles, mat, atoms, heavy, bond, rads, parameters = {}):
    frags = {}
    valence = {'C':4, 'O':2, 'N':3}
    for i, row in enumerate(mat):
        for j, col in enumerate(row):
            if atoms[i] in heavy and atoms[j] in heavy:
               if float(col) > 0:
                   vali = valence[atoms[i]]
                   valj = valence[atoms[j]]
                   if i in rads:
                       vali -= 1
                   if j in rads:
                       valj -= 1
                   vali -= int(round(float(col)))
                   valj -= int(round(float(col)))
                   fragi = atoms[i]
                   fragj = atoms[j]
                   saturate = ''
                   if vali == 1:
                       saturate = 'H'
                   elif vali > 1:
                       saturate = 'H{:g}'.format(vali)
                   fragi = '[{}{}]'.format(fragi, saturate)
                   saturate = ''
                   if valj == 1:
                       saturate = 'H'
                   elif valj > 1:
                       saturate = 'H{:g}'.format(valj)
                   fragj = '[{}{}]'.format(fragj, saturate)
                   frag = fragi + bond[str(int(round(float(col))))] + fragj
                   frag = ob.get_slabel(frag).split('_')[0]
                   if frag in frags:
                       frags[frag] += 1
                   else: 
                       frags[frag] = 1
    #BALANCE
    newfrags = {}
    zedfrags = CBHzed(smiles, mat, atoms, heavy, bond, rads)
    new = {}
    frags =  {k: v for k, v in frags.items() if v}
    for frag in frags:
        newfrags[frag] = frags[frag]
        atoms, mat, rads = get_x2zparams(frag, parameters)
        new = CBHzed(frag, mat, atoms, heavy, bond, rads)
        for n in new:
           if n in newfrags:
               newfrags[n] -= new[n] * frags[frag]
           else:
               newfrags[n] =- new[n] * frags[frag]
    if not frags:
        frags = CBHzed(smiles, mat, atoms, heavy, bond, rads)
    for frag in zedfrags:
        if frag in newfrags:
            newfrags[frag] += zedfrags[frag]
        else:
            newfrags[frag]  = zedfrags[frag]
    frags = newfrags
    frags =  {k: v for k, v in frags.items() if v}
    frags = make_balanced(smiles, frags)
    return frags

def CBHtwo(smiles, mat, atoms, heavy, bond, rads):
    frags = {}
    fmat = square_mat(mat)
    valence = {'C':4, 'O':2, 'N':3}
    for i, row in enumerate(fmat):
        if atoms[i] in heavy:
            
            vali = valence[atoms[i]]
            fragi = atoms[i]
            if i in rads:
                vali -= 1
            vali -= bond_sum(fmat, atoms, heavy, i)
            saturate = ''
            if vali == 1:
                saturate = 'H'
            elif vali > 1:
                saturate = 'H{:g}'.format(vali)
            fragi = '[{}{}]'.format(fragi, saturate)
            saturate = ''
            count = 0
            frag = fragi
            for j, col in enumerate(row):
                if atoms[j] in heavy:
                    if float(col) > 0:
                        count += 1
                        valj = valence[atoms[j]]
                        if j in rads:
                            valj -= 1
                        valj -= int(round(float(col)-0.5))
                        fragj = atoms[j]
                        saturate = ''
                        if valj == 1:
                            saturate = 'H'
                        elif valj > 1:
                            saturate = 'H{:g}'.format(valj)
                        fragj = '[{}{}]'.format(fragj, saturate)
                        frag  += '(' + bond[col] + fragj + ')'
                        vali = valence[atoms[i]]
                        fragi = atoms[i]
                        if i in rads:
                            vali -= 1
                        vali -= int(round(float(col)-0.5))
                        saturate = ''
                        if vali == 1:
                            saturate = 'H'
                        elif vali > 1:
                            saturate = 'H{:g}'.format(vali)
                        fragi = '[{}{}]'.format(fragi, saturate)
                        subfrag = fragi +  bond[col] + fragj 
                        if i < j:
                            subfrag = ob.get_slabel(subfrag).split('_')[0]
                            if subfrag in frags:
                                frags[subfrag] -= 1
                            else: 
                                frags[subfrag] = -1
            frag = ob.get_slabel(frag).split('_')[0]
            if frag in frags:
                frags[frag] += 1
            else: 
                frags[frag] = 1
    frags = make_balanced(smiles, frags)
    return frags

def get_conns(i, axis, mat, atoms, heavy, bond):
    fmat = square_mat(mat)
    ret = atoms[i]
    row = fmat[i]
    for j, col in enumerate(row):
        if j != axis and atoms[j] in heavy and float(col) > 0:
             ret += '({0}{1})'.format(bond[col], atoms[j])
    return ret

def square_mat(mat):
    import copy
    fmat = copy.deepcopy(mat)
    for i, row in enumerate(fmat):
        for j, col in enumerate(fmat): 
            if i == j:
                fmat[i].append('0')
            elif i < j:
                fmat[i].append(mat[j][i])
    return fmat    

def get_recentxyz(s, results):
    xyz = ''
    slabel = ob.get_slabel(s)
    if slabel in results:
        for key in results:
            if 'xyz' in results:
                xyz = results[slabel]['xyz']
    return xyz

def choose_xyz(smiles, parameters):
    xyz = ''
    if 'all results' in parameters:
        xyz = get_recentxyz(smiles, parameters['all results'])
    if not xyz:
        if 'xyzdir' in parameters:
           xyz = io.read_xyzdir(smiles, parameters['xyzdir'])
    if not xyz:
            xyz = ob.get_xyz(smiles)
    return xyz

def get_x2zparams(smiles, parameters={}):
    """
    returns relevant x2z information for CBH calculation
    INPUT: 
    parameters -- dicitonary with qtc data (can be empty)
    OUTPUT:
    atoms      -- list of atoms in order of the x2z adjacency matrix
    mat        -- triangular connectivity matrix from x2z
    rads       -- list of radical sites zero indexed that correlate to the order of the atoms list
    EXAMPLE:
    >> get_x2zparams('C[C]=O',{})
    ['C', 'C', 'O', 'H', 'H', 'H'], [[], ['1'], ['0', '2'], ['1', '0', '0'], ['1', '0', '0', '0'], ['1', '0', '0', '0', '0']],  [1]
    """
    xyz = choose_xyz(smiles, parameters)
    io.write_file(xyz,'x2z.xyz')
    lines = qc.run_x2z('x2z.xyz', 'x2z')
    atoms, mat, rads = get_adj_mat(lines)
    return atoms, mat, rads

def cbh_coefficients(smiles, ref, parameters = {}):
   
    logging.info(smiles)
    if ob.get_natom(smiles) < 2:
        clist = [1]
        fraglist = [smiles]
        msg = 1
    else:
        atoms, mat, rads = get_x2zparams(smiles, parameters)
        heavy = ['C','O','N']
        bond  = {'1':'','2':'=','3':'#','4':'$','1.5':':','2.5':'=','1.7':'=','2.3':'=','1.3':'','2.7':'#'}
        
        msg = ''
        if ref.lower() == 'cbh0':
            frags = CBHzed(smiles, mat, atoms, heavy, bond, rads)
        elif ref.lower() == 'cbh1':
            frags = CBHone(smiles, mat, atoms, heavy, bond, rads, parameters)
        elif ref.lower() == 'cbh2':
            frags = CBHtwo(smiles, mat, atoms, heavy, bond, rads)

        output = print_lhs_rhs(smiles, frags)
        if isUnbalanced(get_balance(smiles, frags)):
            output += '\nERROR: unbalanced fragmentation'
        msg += output
        msg += '\n'

        fraglist = []
        clist = []
        for frag in frags:
            fraglist.append(frag)
            clist.append(frags[frag])
    return clist, fraglist, msg

def get_bac(parameters, mylist, samppercent = 0, errthresh = 10):
    """
    Computes BAC parameters (corrections to the HoF for each bond type for each level of theory)
    INPUT
    parameters  -- dictionary containing species info
    mylist      -- list of species
    samppercent -- percentage of species you want to form the set with (0 being all, 1 being none)
    errthresh   -- threshold for what species will be ignored when forming BAC based on kj difference from ANL0 
    OUTPUT
    out         -- msg to log
    bac.txt     -- BAC parameters file
    """
    out = 'begin least squares\n'
    bondlabel = []
    mols = []
    taskmat = {}
    testset = ''
    bacset = ''
    for i,s in enumerate(mylist):
        anl = nest_2_dic(qc.get_slabel(s), 'delHf','ANL0')
        bonds = None
        if anl:
            if np.random.rand() > samppercent:
                sresults = parameters['all results'][qc.get_slabel(s)]
                for qcresultkey, qcresultval in sorted(sresults.iteritems(),key= lambda x: x[0]):
                    logging.debug(s)
                    task = qcresultkey#qcresultval['chemkin'].splitlines()[0].split()[2]
                    if task in taskmat:
                        mat = taskmat[task]['mat']
                        err = taskmat[task]['errors']
                    else:
                        mat = []
                        err = []
                        taskmat[task] = {}
                        taskmat[task]['mat'] = mat
                        taskmat[task]['errors'] = err
                    energy = qcresultval['deltaH0']
                    error = anl / ut.kcal2kj - energy
                    #error = energy - anl/ ut.kcal2kj
                    if abs(error) < errthresh:
                        mols.append(qc.get_slabel(s))
                        natom = ob.get_natom(ob.get_mol(s))
                        if not bonds:
                            if natom > 1:
                                xyz = ''
                                if 'xyz' in qcresultval:
                                    xyz = qcresultval['xyz']
                                if not xyz:
                                    xyz = ob.get_xyz(ob.get_mol(s))
                                x2zout = qc.run_x2z(xyz, parameters['x2z'])
                                bacset += '{}\n'.format( qc.get_slabel(s)  )
                                x2zout = qc.run_x2z(xyz, parameters['x2z'])
                                bonds =  qc.get_x2z_bonds(x2zout)
                            else:
                                bonds = {}
                            hofbasis = qcresultval['heat of formation basis']
                            hofcoeff = qcresultval['heat of formation coeff']
                            hofbonds = {}
                            for c, bas in enumerate(hofbasis):
                                nnatom = ob.get_natom(ob.get_mol(bas))
                                if nnatom > 1:
                                    xyz = ''
                                    if 'xyz' in parameters['all results'][qc.get_slabel(bas)][task]:
                                        xyz = parameters['all results'][qc.get_slabel(bas)][task]['xyz']
                                    if not xyz:
                                        xyz = ob.get_xyz(ob.get_mol(s))
                                    x2zout = qc.run_x2z(xyz, parameters['x2z'])
                                    hofbonds =  qc.get_x2z_bonds(x2zout)
                                    for bond in hofbonds:
                                        if bond in bonds:
                                            bonds[bond] -= hofcoeff[c] * hofbonds[bond]
                                        else:
                                            bonds[bond] =- hofcoeff[c] * hofbonds[bond]
                        bondarray = np.zeros(100)
                        for bond in bonds:
                            new = True
                            j = 0
                            for label in bondlabel:
                                if bond == label:
                                    bondarray[j] = bonds[bond]
                                    new = False 
                                    break
                                j += 1
                            if new:
                                bondlabel.append(bond)
                                bondarray[j] = bonds[bond]
                        mat.append(bondarray.tolist())
                        err.append(error)
                        taskmat[task]['mat'] = mat
                    else:
                        out += 'ommitting {} from BAC, error is enormous\n'.format(s)
                    taskmat[task]['errors'] = err
            else:
                testset += '{}\n'.format(qc.get_slabel(s))
    io.write_file(bacset, 'formset.txt')
    io.write_file(testset, 'testset.txt')
    bacout = 'Formed from: {}\n'.format( ' ,'.join(mols)) 
    for task in taskmat:
        mat = taskmat[task]['mat']
        mat = [mat[k][:len(bondlabel)] for k in range(len(mat))]
        errors =  taskmat[task]['errors']
        findbac = True
        try:
            #print task
            #print 'AVG \t {:.4f} \nRMS \t {:.4f}'.format(np.average(errors),np.sqrt(np.average(np.square(errors))))
            coeff = np.linalg.lstsq(mat, errors)[0]
            bacout +=  '\n' + task + '\n'
            for k, label in enumerate(bondlabel):
                bacout += '{}  {:4.3f}\n'.format(label, coeff[k])
            out += 'BAC parameters successfully computed for {}, appended to bac.txt\n'.format(task)
        except:
           out += 'Cannot compute BAC for {}, nan found in errors\n'.format(task)
    io.write_file(bacout, 'strictestbac.txt')
    return out
      
def calc_bac(parameters, bonds, task):
    """
    Given a set of bonds, calculate the bond additivity correction
    using a  template specified in parameters
    """
    correction = 0.
    bacfile = parameters['bacdirectory'] +'/bac.txt'
    #bacfile = 'bac.txt'
    bacfile = io.read_file(bacfile)
    if task + '\n' in bacfile:
        bacfile = bacfile.split(task + '\n')[1].split('\n\n')[0]
        bacfile = bacfile.splitlines()
        bacdic  = {}
        for line in bacfile:
            line = line.split()
            if len(line) > 1:
                bacdic[line[0]] = float(line[1])
        for key in bonds:
            if key in bacdic:
                correction += bonds[key] * bacdic[key] 
            else:
                logging.info( '{} bond has no correction'.format(key))
    return correction



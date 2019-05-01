def read_freqs(file_str):
    """ freqs
    """
    freqs = []
    for line in file_str.splitlines():
        if 'i' not in line:
            freqs.append(float(line.strip()))                
        
    return freqs


def read_sym_factor(file_str):
    """ symmetry factor
    """
    sym_factor = float(file_str.strip())
    return sym_factor
   

def read_elec_levels(file_str):
    """ elec_levels
    """
    elec_levels = []
    for line in file_str.splitlines():
        elec_levels.append([float(val) for val in line.strip().split()])
        
    return elec_levels


def read_hindered_rotors(file_str):
    """ hindered rotors; calls from several files
    """

    hindlines = file_str.splitlines()
    hind_rot = []
    for i in range(len(hindlines)):
        if 'Rotor' in hindlines[i] and 'Hindered' in hindlines[i]:
            group =     [ int(val) for val in hindlines[i+1].split()[1:] ]
            axis =      [ int(val) for val in hindlines[i+2].split()[1:] ]
            symmetry =  int(hindlines[i+3].split()[1])
            potential = [ float(val) for val in hindlines[i+5].split() ]
            hind_rot.append( [group, axis, symmetry, potential] )            
    
    return hind_rot

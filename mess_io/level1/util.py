def read_freqs(file_str):
    """ freqs
    """
    freqs = []
    for line in file_str.splitlines():
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



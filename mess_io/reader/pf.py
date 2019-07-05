"""
  Reads the output of messpf
"""

def read_pf(pf_file_name):
    """ reads the pf.dat file.
    """
    
    temps, logQ, dQ, dQ2 = [], [], [], []
    with open(pf_file_name, 'r') as pf_file:
        for i, line in enumerate(pf_file):
            if i != 0:
                tmp = line.strip().split()
                temps.append(tmp[0]) 
                logQ.append(tmp[1]) 
                dQ.append(tmp[2]) 
                dQ2.append(tmp[3]) 

    return temps, logQ, dQ, dQ2

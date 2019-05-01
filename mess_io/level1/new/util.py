"""
libraries of functions to write the mess input
"""

import os
import autofile
import mess_io.writer


# FUNCTIONS TO HELP BUILD MORE COMPLICATED PARTS OF THE MESS INPUT #

def build_core(core_type,
               geom1='', geom2='', stoich='',
               sym_factor='', 
               ne_file='',
               interp_emax=100, quant_lvl_emax=9,
               pot_prefactor=10, pot_power=6):
    """  writes the core section since it is a bit more complicated 
         might want to break up later, but this works for now
    """
    
    # Get the string for the core using the geometry
    if core_type == 'rigidrotor':
        core_str = mess_io.writer.write_core_rigidrotor(geom1, sym_factor)
    elif core_type == 'multirotor':
        core_str = mess_io.writer.write_core_multirotor(geom1, sym_factor, pot_surf, rotor_int_str,
                                                        interp_emax=interp_emax, 
                                                        quant_lvl_emax=quant_lvl_emax)
    elif core_type == 'phasespace':
        core_str = mess_io.writer.write_core_phasespace(geom1, geom2, sym_factor, stoich,
                                                        pot_prefactor=pot_prefactor, 
                                                        pot_power_exp=pot_power_exp)
    elif core_type == 'rotd':
        core_str = mess_io.writer.write_core_rotd(sym_factor, ne_file, stoich)
    else:
        raise NotImplementedError

    return core_str


def build_hr(hr_list):
    """ builds a MESS input string for a sequence of hindered rotors 
    """

    hr_str = ''
    for hr in hr_list:
        group, axis, symmetry, potential = hr[0], hr[1], hr[2], hr[3]
        hr_str += mess_io.writer.write_rotor_hindered(group, axis, symmetry, potential)

    return hr_str


# FUNCTIONS TO GET DATA FROM DIRECTORY PATHS #

def energy_from_path(ref_elec=('', ''), ref_zpve=('', ''),
                     spec1_elec=('', ''), spec1_zpve=('', ''),
                     spec2_elec=('', ''), spec2_zpve=('', '')):
    """ obtains a relative energy for a unimolecular or bimolecular species
        from a series of paths
    """

    # Read in the reference electronic energy and zpve
    ref_elec_path = os.path.join(ref_elec[0], ref_elec[1])
    ref_zpve_path = os.path.join(ref_zpve[0], ref_zpve[1])
    with open(ref_elec_path, 'r') as f:
        ref_energy_str = f.read()
    with open(ref_zpve_path, 'r') as f:
        ref_zpve_str = f.read()
    ref_e_elec = autofile.read.energy(ref_energy_str)
    ref_e_zpve = autofile.read.energy(ref_zpve_str)

    # Read in the species 1 electronic energy and zpve
    spec1_elec_path = os.path.join(spec1_elec[0], spec1_elec[1])
    spec1_zpve_path = os.path.join(spec1_zpve[0], spec1_zpve[1])
    with open(spec1_elec_path, 'r') as f:
        energy_str1 = f.read()
    with open(spec1_zpve_path, 'r') as f:
        zpve_str1 = f.read()
    e_elec1 = autofile.read.energy(energy_str1)
    e_zpve1 = autofile.read.energy(zpve_str1)
    
    # Read in the species 2 electronic energy and zpve
    # Values will be set to 0.0 if rel energy is being computed for a single species
    if spec2_elec == ('', ''):
        e_elec2 = 0.0
    else: 
        spec2_elec_path = os.path.join(spec2_elec[0], spec2_elec[1])
        with open(spec2_elec_path, 'r') as f:
            energy_str2 = f.read()
        e_elec2 = autofile.read.energy(energy_str2)
    if spec2_zpve == ('', ''):
        e_zpve2 = 0.0
    else: 
        spec2_zpve_path = os.path.join(spec2_zpve[0], spec2_zpve[1])
        with open(spec2_zpve_path, 'r') as f:
            zpve_str2 = f.read()
        e_zpve2 = autofile.read.energy(zpve_str2)

    # Calculate the reference and species energies
    ref_ene = ref_e_elec + ref_e_zpve
    ene1 = e_elec1 + e_zpve1
    ene2 = e_elec2 + e_zpve2

    # Compute the rel energy
    rel_energy = ((ene1 + ene2) - ref_ene) * 627.5095

    return rel_energy


def geom_from_path(file_path, file_name):
    """ obtains a geometry from a path
    """

    # Set the file path and read in the file string    
    geom_file = os.path.join(file_path, file_name)
    with open(geom_file, 'r') as f:
        geom_str = f.read()
    
    # Obtain a geom object from the string
    geom = autofile.read.geometry(geom_str)

    return geom


def freqs_from_path(file_path, file_name):
    """ obtains a frequencies from a path
    """
    
    # Set the file path and read in the file string    
    freqs_file = os.path.join(file_path, file_name)
    with open(freqs_file, 'r') as f:
        freqs_str = f.read()
    
    # Obtain a freqs object from the string
    freqs = read_freqs(freqs_str)

    return freqs


def hr_from_path(file_path, file_name):    
    """ obtains hindered rotors from a path
    """

    # Set the file path and read in the file string    
    hr_file = os.path.join(file_path, file_name)
    with open(hr_file, 'r') as f:
        hr_str = f.read()
    
    # Obtain a hindered-rotors object from the string
    hrs = read_hindered_rotors(hr_str)

    return hrs


def eckart_from_path(ts_freqs_path=('', ''), 
                     ts_energy_path=('', ''),
                     reac_energy_path=('', ''), 
                     prod_energy_path=('', '')):
    """ Get Eckart tunneling
    """

    # Get the imaginary frequency
    ts_freqs_file = os.path.join(ts_freqs_path[0], ts_freqs_path[1])
    with open(ts_freq_file, 'r') as f:
        ts_freqs_str = f.read()
    imag_freq = read_imag_freq(ts_freqs_str)

    # Get the energies to compute well depths
    ts_file = os.path.join(ts_path[0], ts_path[1])
    reac_file = os.path.join(reac_path[0], reac_path[1])
    prod_file = os.path.join(prod_path[0], prod_path[1])
    with open(ts_file, 'r') as f:
        ts_str = f.read()
    e_ts = autofile.read.energy(ts_str)
    with open(reac_file, 'r') as f:
        reac_str = f.read()
    e_reac = autofile.read.energy(reac_str)
    with open(prod_file, 'r') as f:
        prod_str = f.read()
    e_prod = autofile.read.energy(prod_str)
   
    # Calculate the well depths
    well_depth1 = ( e_ts - e_reac ) * 627.5095         
    well_depth2 = ( e_ts - e_prod ) * 627.5095         
    e_elec1 = autofile.read.energy(energy_str1)


    return imag_freq, well_depth1, well_depth2


# FUNCTIONS TO READ DATA FROM FILES, COUPLED TO ABOVE PATH FUNCTIONS; SHOULD BE REPLACED #

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
    """ hindered rotors
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


def read_imag_freq(file_str):
    """ imag freq
    """
    for line in file_str.splitlines():
        if 'i' in line:
            imag_freq = line.strip().replace('i', '')
        
    return imag_freq

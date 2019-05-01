"""
libraries of functions to write the mess input
"""

import os
import autofile
import mess_io.writer
import util


# FOR SECTIONS #

def global_keys(filename='', messtype='', temperatures='', pressures=''):
    """ Writes the global keys section
    """
    
    # Use the writer to create a string for the global keyword section
    if messtype == 'pf':
        global_keys_str = mess_io.writer.write_global_pf(
                            temperatures, rel_temp_inc=0.001, atom_dist_min=0.6)
    elif messtype == 'reaction':
        global_keys_str = mess_io.writer.write_global_reaction(
                            temperatures, pressures)
    else:
        raise NotImplementedError

    # Write the global section string
    with open(filename, 'a') as f:
        f.write(global_keys_str)


def energy_transfer(filename='',
                    exp_factor='', exp_power='', exp_cutoff='',
                    eps1='', eps2='',
                    sig1='', sig2='',
                    mass1='', mass2=''):
    """ Write the energy transfer section
    """

    # Use the writer to create a string for the energy transfer section
    energy_trans_str = mess_io.writer.write_energy_transfer(exp_factor, exp_power, exp_cutoff,
                                                            eps1, eps2,
                                                            sig1, sig2,
                                                            mass1, mass2)
    
    # Write the energy transfer section string
    with open(filename, 'a') as f:
        f.write(energy_trans_str)


def species(filename='', label='', data=''):
    """ species writer
    """

    # Write the species string
    species_str = mess_io.writer.write_species(label, data)

    # Write the species section string
    with open(filename, 'a') as f:
        f.write(species_str)


def well(filename='', label='', data=''):
    """ well writer
    """

    # Write the species string
    well_str = mess_io.writer.write_well(label, data)

    # Write the well section string
    with open(filename, 'a') as f:
        f.write(well_str)


def bimolecular(filename='', bimol_label='', 
                species1_label='', species1_data='',
                species2_label='', species2_data='',
                ground_energy=''):
    """ bimolecular writer
    """

    # Write the bimolecular string
    bimolecular_str = mess_io.writer.write_bimolecular(bimol_label,
                                                       species1_label, species1_data,
                                                       species2_label, species2_data,
                                                       ground_energy)

    # Write the bimolecular section string
    with open(filename, 'a') as f:
        f.write(bimolecular_str)


def ts_sadpt(filename='',
             ts_label='', reac_label='', prod_label='',
             ts_data=''):
    """ ts saddle-point writer
    """

    # Write the ts saddle-point string
    ts_sadpt_str = mess_io.writer.write_ts_sadpt(ts_label, reac_label, prod_label, ts_data)
    
    # Write the ts saddle-point section string
    with open(filename, 'a') as f:
        f.write(ts_sadpt_str)


# SECTION HEADER and SEPERATOR STRINGS #

def rxn_chan_head(filename=''):
    """ writes a header for the reaction channel section
    """
    
    # Get the header string from the library
    head_str = mess_io.writer.stringslib.RXN_CHAN_HEAD_STR

    # Write the header string
    with open(filename, 'a') as f:
        f.write(head_str)


def species_head(filename=''):
    """ writes a header for the species section
    """
    
    # Get the header string from the library
    head_str = mess_io.writer.stringslib.SPECIES_HEAD_STR

    # Write the header string
    with open(filename, 'a') as f:
        f.write(head_str)


def species_sep(filename=''):
    """ writes a string for seperating species sections
    """
    
    # Get the header string from the library
    sep_str = mess_io.writer.stringslib.SPECIES_SEC_SEP_STR

    # Write the header string
    with open(filename, 'a') as f:
        f.write(sep_str)


# MOLECULE AND ATOM SECTION WRITERS #

def atom(name='', elec_levels=''): 
    """ atom writer
    """

    # Use the writer to create a string for the atom section
    atom_str = mess_io.writer.write_atom(name, elec_levels)

    return atom_str


def molecule(core='', freqs='', zero_energy='', elec_levels='',
             hind_rot='', tunnel='',
             anharm='', rovib_coups='', rot_dists=''):
    """ molecule writer
    """

    # Use the writer to create a string for the molecule section
    # can't do sct_tunneling, anharm, rovib_coups, and rot_dists reading
    molecule_str = mess_io.writer.write_molecule(
        core, freqs, zero_energy, elec_levels,
        hind_rot=hind_rot, tunnel=tunnel,
        anharm=anharm, rovib_coups=rovib_coups, rot_dists=rot_dists)

    return molecule_str


def core(core_type,
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


def hr(hr_list):
    """ tries to 
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
    if spec2_elec == ('', ''):
        e_elec2 = 0.0
    else: 
        spec1_elec_path = os.path.join(spec1_elec[0], spec1_elec[1])
        with open(spec2_elec_path, 'r') as f:
            energy_str2 = f.read()
        e_elec2 = autofile.read.energy(energy_str2)
    if spec2_zpve == ('', ''):
        e_zpve2 = 0.0
    else: 
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
    
    geom_file = os.path.join(file_path, file_name)
    with open(geom_file, 'r') as f:
        geom_str = f.read()
    geom = autofile.read.geometry(geom_str)

    return geom


def freqs_from_path(file_path, file_name):
    """ obtains a frequencies from a path
    """
    
    freq_file = os.path.join(file_path, file_name)
    with open(freqs_path, 'r') as f:
        freqs_str = f.read()
    freqs = util.read_freqs(freqs_str)

    return freqs


def hr_from_path(file_path, file_name):    
    """ obtains hindered rotors from a path
    """

    hr_file = os.path.join(file_path, file_name)
    with open(hr_file, 'r') as f:
        hr_str = f.read()
    hrs = util.read_hindered_rotors(hr_str)

    return hrs

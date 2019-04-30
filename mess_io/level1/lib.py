"""
libraries of functions to write the mess input
"""

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


def molecule(core='', 
             zero_energy='',
             geom='', 
             sym_factor='',
             freqs='',
             elec_levels=''):
    """ molecule writer
    """

    # Get the string for the core using the geometry
    if core == 'rigidrotor':
        core_str = mess_io.writer.write_core_rigidrotor(geom, sym_factor)
    
    # Use the writer to create a string for the molecule section
    molecule_str = mess_io.writer.write_molecule(
        core_str, freqs, zero_energy, elec_levels) 

    return molecule_str


# FUNCTIONS TO GET DATA FROM DIRECTORY PATHS #

def energy_from_path(ref_elec_path='', ref_zpve_path='',
                     spec1_elec_path='', spec1_zpve_path='',
                     spec2_elec_path='', spec2_zpve_path=''):
    """ obtains a relative energy for a unimolecular or bimolecular species
        from a series of paths
    """
   
    # Read in the reference electronic energy and zpve
    with open(ref_elec_path, 'r') as f:
        ref_energy_str = f.read()
    with open(ref_zpve_path, 'r') as f:
        ref_zpve_str = f.read()
    ref_e_elec = autofile.read.energy(ref_energy_str)
    ref_e_zpve = autofile.read.energy(ref_zpve_str)

    # Read in the species 1 electronic energy and zpve
    with open(spec1_elec_path, 'r') as f:
        energy_str1 = f.read()
    with open(spec1_zpve_path, 'r') as f:
        zpve_str1 = f.read()
    e_elec1 = autofile.read.energy(energy_str1)
    e_zpve1 = autofile.read.energy(zpve_str1)
    
    # Read in the species 2 electronic energy and zpve
    if spec2_elec_path == '':
        e_elec2 = 0.0
    else: 
        with open(spec2_elec_path, 'r') as f:
            energy_str2 = f.read()
        e_elec2 = autofile.read.energy(energy_str2)
    if spec2_zpve_path == '':
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


def geom_from_path(geom_path):
    """ obtains a geometry from a path
    """
    
    with open(geom_path, 'r') as f:
        geom_str = f.read()
    geom = autofile.read.geometry(geom_str)

    return geom


def freqs_from_path(freqs_path):
    """ obtains a frequencies from a path
    """
    
    with open(freqs_path, 'r') as f:
        freqs_str = f.read()
    freqs = util.read_freqs(freqs_str)

    return freqs

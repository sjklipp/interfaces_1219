"""
libraries of functions to write the mess input
"""

import autofile
import mess_io.writer
import util

# FOR SECTIONS #

def globalkeys(filename='', messtype='', temperatures='', pressures=''):
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


def bimolecular(filename='', label='', 
                species1_label='', species1_data='',
                species2_label='', species2_data=''
                ref_energy_path=''):
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



# SECTION HEADER STRINGS #

def species_head(filename=''):
    """ writes a header for the species head
    """
    
    # Get the header string from the library
    head_str = mess_io.writer.stringslib.SPECIES_HEAD_STR

    # Write the header string
    with open(filename, 'a') as f:
        f.write(head_str)


# MOLECULE AND ATOM DATA #

def atom_from_data(name='', elec_levels=''): 
    """ atom writer
    """

    # Use the writer to create a string for the atom section
    atom_str = mess_io.writer.write_atom(name, elec_levels)

    return atom_str


def atom_from_path(name='', elec_levels=''): 
    """ atom writer
    """
    
    # Elec Levels:
    with open(elec_levels_path, 'r') as f:
        elec_levels_str = f.read()
    elec_levels = util.read_elec_levels(elec_levels_str)

    # Use the writer to create a string for the atom section
    atom_str = mess_io.writer.write_atom(name, elec_levels)

    return atom_str


def molecule_from_data(core='', 
                       geom='', 
                       symfactor='', 
                       freqs='', 
                       elec_levels='', 
                       zero_e=''):
    """ molecule writer
    """

    # Get the string for the core using the geometry
    if core == 'rigidrotor':
        core_str = mess_io.writer.write_core_rigidrotor(geom, symfactor)
    
    # Use the writer to create a string for the molecule section
    molecule_str = mess_io.writer.write_molecule(
        core_str, freqs, zero_e, elec_levels) 

    return molecule_str


def molecule_from_path(core='', 
                       geom_path='', 
                       energy_path='',
                       ref_energy_path='',
                       freqs_path='',
                       sym_factor_path='',
                       elec_levels_path=''):
    """ molecule writer
    """

    # Get the pieces of information from the path
    # Geometry:
    with open(geom_path, 'r') as f:
        geom_str = f.read()
    geom = autofile.read.geometry(geom_str)
    
    # Energy:
    with open(energy_path, 'r') as f:
        energy_str = f.read()
    ene = autofile.read.energy(energy_str)
    with open(ref_energy_path, 'r') as f:
        ref_energy_str = f.read()
    ref_ene = autofile.read.energy(ref_energy_str)
    zero_e = (ene - ref_ene) * 627.5095
    
    # Symmetry Factor:
    with open(sym_factor_path, 'r') as f:
        sym_factor_str = f.read()
    sym_factor = util.read_sym_factor(sym_factor_str)

    # Frequencies:
    with open(freqs_path, 'r') as f:
        freqs_str = f.read()
    freqs = util.read_freqs(freqs_str)

    # Elec Levels:
    with open(elec_levels_path, 'r') as f:
        elec_levels_str = f.read()
    elec_levels = util.read_elec_levels(elec_levels_str)


    # Get the string for the core using the geometry
    if core == 'rigidrotor':
        core_str = mess_io.writer.write_core_rigidrotor(geom, sym_factor)
    
    # Use the writer to create a string for the molecule section
    molecule_str = mess_io.writer.write_molecule(
        core_str, freqs, zero_e, elec_levels) 

    return molecule_str

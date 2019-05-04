"""
mess molecule string
"""

import mess_io
from mess_io.writer_lvl1.params import MESS
from mess_io.writer_lvl1.util import geom_from_path
from mess_io.writer_lvl1.util import energy_from_path
from mess_io.writer_lvl1.util import freqs_from_path
from mess_io.writer_lvl1.util import hr_from_path
from mess_io.writer_lvl1.util import build_hr


def build_mess_molecule_str(data_dict, path_dict):
    """ builds the string using the species dictionaries
    """

    # Geometry set directly from dictionary or using paths
    if data_dict[MESS.GEOM] == 'PATH':
        geom = geom_from_path(path_dict[MESS.MOL_PATH], path_dict[MESS.GEOM_PATH])    
    else:
        geom = data_dict[MESS.GEOM]
    
    # Symmetry Factor set directly from dictionary or using paths
    if data_dict[MESS.SYMF] == 'PATH':
        sym_factor = sym_factor_from_path(path_dict[MESS.MOL_PATH], path_dict[MESS.SYMF_PATH])    
    else:
        sym_factor = data_dict[MESS.SYMF]
    
    # Frequencies set directly from dictionary or using paths
    if data_dict[MESS.FREQ] == 'PATH':
        freqs = freqs_from_path(path_dict[MESS.MOL_PATH], path_dict[MESS.FREQ_PATH])    
    else:
        freqs = data_dict[MESS.FREQ]
    
    # Electronic Levels set directly from dictionary or using paths
    if data_dict[MESS.ELVL] == 'PATH':
        elec_levels = elec_levels_from_path(path_dict[MESS.MOL_PATH], path_dict[MESS.ELVL_PATH])    
    else:
        elec_levels = data_dict[MESS.ELVL]
   
    # Zero-Energy set directly from dictionary or using paths; species2 sections of fxn call set to ''
    if data_dict[MESS.ZENE] == 'PATH':
        # Get the total energy for the reference    
        ref_total_energy = energy_from_path(path_dict[MESS.REF_PATH], 
                                            path_dict[MESS.REF_ENE_PATH], 
                                            path_dict[MESS.REF_ZPVE_PATH])
        # Get the total energy for the species 1    
        spec_total_energy = energy_from_path(path_dict[MESS.MOL_PATH], 
                                             path_dict[MESS.MOL_ENE_PATH], 
                                             path_dict[MESS.MOL_ZPVE_PATH])
        # Calculate the relative ground energy using the species energies
        zero_energy = spec_total_energy - ref_total_energy
        zero_energy *= 627.5095
    else:
        zero_energy = data_dict[MESS.ZENE]

    # SET ADDITIONAL KEYWORDS (AVAIL BY PATH OR VALUE) IF THEY ARE REQUESTED
    
    # Hindered Rotors set directly from dictionary or using paths
    if MESS.HROT in data_dict:
        if data_dict[MESS.HROT] == 'PATH':
            hind_rot_str = build_hr(hr_from_path(path_dict[MESS.MOL_PATH], path_dict[MESS.HROT_PATH]))    
        else:
            hind_rot_str = build_hr([data_dict[MESS.HROT]])
    else:
            hind_rot_str = ''
    
    # Anharmonicity set directly from dictionary or using paths
    if MESS.ANHM in data_dict:
        if data_dict[MESS.ANHM] == 'PATH':
            anharm = anharm_from_path()
        else:
            anharm = data_dict[MESS.ANHM]
    else:
            anharm = ''
    
    # Rovibrational Coupling set directly from dictionary or using paths
    if MESS.RVCP in data_dict:
        if data_dict[MESS.RVCP] == 'PATH':
            rovib_coups = rovib_coups_from_path()
        else:
            rovib_coups = data_dict[MESS.RVCP]
    else:
            rovib_coups = ''
    
    # Rotational Distortion set directly from dictionary or using paths
    if MESS.RDIS in data_dict:
        if data_dict[MESS.RDIS] == 'PATH':
            rot_dists = rot_dist_from_path()
        else:
            rot_dists = data_dict[MESS.RDIS]
    else:
            rot_dists = ''
    
    # Tunneling set directly from dictionary or using paths
    if MESS.TUNL in data_dict:
        if data_dict[MESS.TUNL][0] == 'ECKART' and data_dict[MESS.TUNL][1] == 'PATH':
            
            # Get the imaginary frequency
            imag_freq = imag_freq_from_path(path_dict[MESS.MOL_PATH], path_dict[MESS.FREQ_PATH])
            
            # Get the total energy for the reference    
            ts_total_energy = energy_from_path(path_dict[MESS.REF_PATH], 
                                               path_dict[MESS.REF_ENE_PATH], 
                                               path_dict[MESS.REF_ZPVE_PATH])
            
            # Get the total energy for the reactant
            reac_total_energy = energy_from_path(path_dict[MESS.REAC_PATH], 
                                                 path_dict[MESS.REAC_ENE_PATH], 
                                                 path_dict[MESS.REAC_ZPVE_PATH])
            
            # Get the total energy for the product
            prod_total_energy = energy_from_path(path_dict[MESS.PROD_PATH], 
                                                 path_dict[MESS.PROD_ENE_PATH], 
                                                 path_dict[MESS.PROD_ZPVE_PATH])
            
            # Get the reactant and product well depths
            well_depth1 = ts_total_energy - reac_total_energy
            well_depth1 *= 627.5096
            well_depth2 = ts_total_energy - prod_total_energy
            well_depth2 *= 627.5096

            # Get the imaginary frequency
            imag_freq = imag_freq_from_path(path_dict[MESS.MOL_PATH], path_dict[MESS.FREQ_PATH])

            # Write the eckart tunneling string
            tunnel_str = mess_io.writer.write_tunnel_eckart(imag_freq, well_depth1, well_depth2)
    
        elif data_dict[MESS.TUNL][0] == 'ECKART' and data_dict[MESS.TUNL][1] != 'PATH':

            assert len(data_dict[MESS.TUNL][1]) == 3

            # Get the values needed for Eckart tunneling out of the list
            imag_freq = data_dict[MESS.TUNL][1][0] 
            well_depth1 = data_dict[MESS.TUNL][1][1] 
            well_depth2 = data_dict[MESS.TUNL][1][2] 
            
            # Write the eckart tunneling string
            tunnel_str = mess_io.writer.write_tunnel_eckart(imag_freq, well_depth1, well_depth2)
    else:
            tunnel_str = ''

    # SET ADDITIONAL KEYWORDS (AVAIL BY PATH OR VALUE) IF THEY ARE REQUESTED
    
    # Setting various optional keywords if they are set
    if MESS.INTRP_EMAX in data_dict:
        interp_emax = data_dict[MESS.INTRP_EMAX]
    else:
        interp_emax = 10
    if MESS.QUANT_LVL_EMAX in data_dict:
        quant_lvl_emax = data_dict[MESS.QUANT_LVL_EMAX]
    else:
        quant_lvl_emax = 10
    
    # Build the core string
    if data_dict[MESS.CORE] == 'rigidrotor':
        core_str = mess_io.writer.write_core_rigidrotor(geom, sym_factor)
    #elif data_dict[MESS.CORE] == 'multirotor':
    #    core_str = mess_io.writer.write_core_multirotor(geom, sym_factor, pot_surf, rotor_int_str,
    #                                                    interp_emax=interp_emax, 
    #                                                    quant_lvl_emax=quant_lvl_emax)
    
    # Build the molecule string with the core string and above data
    molecule_str = mess_io.writer.write_molecule(
        core_str, freqs, zero_energy, elec_levels,
        hind_rot=hind_rot_str, tunnel=tunnel_str,
        anharm=anharm, rovib_coups=rovib_coups, rot_dists=rot_dists)
   
    return molecule_str

"""
uses dictionaries defined in other file to build the MESS strings
"""

import os
import mess_io
from mess_io.writer_lvl1.params import MESS
from mess_io.writer_lvl1.molecule_ import build_mess_molecule_str


def build_mess_ts_str(trans_states):
    """ Builds the MESS input strings for each of the species
    """
    
    # Initialize the species string with the species header
    #trans_states_str = mess_io.writer.stringslib.TS_HEAD_STR
    trans_states_str = ''

    # Loop over each transition state
    for trans_state in trans_states:
    
        # Set dictionaries for clarity
        head_dict = trans_state[0]
        data_dict = trans_state[1]
        path_dict = trans_state[2]
    
        # Write string if it is a saddle-point or IRC
        if head_dict[MESS.TS_TYP] == 'sadpt': 
            trans_state_str = build_mess_ts_sadpt_str(head_dict, data_dict, path_dict)
        elif head_dict[MESS.TS_TYP] == 'irc': 
            trans_state_str = build_mess_ts_irc_str(head_dict, data_dict, path_dict)
        else:
            raise NotImplementedError

        # Add to the species string the specie string and a seperator
        trans_states_str += trans_state_str
        trans_states_str += mess_io.writer.stringslib.SPECIES_SEC_SEP_STR
    
    return trans_states_str


def build_mess_ts_sadpt_str(head_dict, data_dict, path_dict):
    """ ts for a saddle point
    """

    # SET THE VALUES FOR REQUIRED INFORMATION #
    
    # Check for existance of required keys
    assert MESS.HD_LBL in head_dict
    assert MESS.R_LBL in head_dict
    assert MESS.P_LBL in head_dict
    assert MESS.GEOM in data_dict
    assert MESS.SYMF in data_dict
    assert MESS.FREQ in data_dict
    assert MESS.ZENE in data_dict
    assert MESS.ELVL in data_dict
    assert MESS.CORE in data_dict
    
    # Check that the core value is appropriate for species
    assert (data_dict[MESS.CORE] == 'rigidrotor' or data_dict[MESS.CORE] == 'multirotor' or
            data_dict[MESS.CORE] == 'rotd' or data_dict[MESS.CORE] == 'phasespace')

    # Build the molecule string with the core string and above data
    molecule_str = build_mess_molecule_str(head_dict, data_dict, path_dict)
    
    # Build the species string with the label and molecule string
    trans_state_str = mess_io.writer.write_ts_sadpt(head_dict[MESS.HD_LBL],
                                                    head_dict[MESS.R_LBL],
                                                    head_dict[MESS.P_LBL],
                                                    molecule_str)

    return trans_state_str


def build_mess_ts_irc_str(head_dict, data_dict, path_dict):
    """ ts for an irc
        only supported:
            rigidrotor core
            symmetrynumber 1.000
    """

    # SET THE VALUES FOR REQUIRED INFORMATION #
    
    # Check for existance of required keys
    assert MESS.HD_LBL in head_dict
    assert MESS.R_LBL in head_dict
    assert MESS.P_LBL in head_dict
    assert MESS.GEOM in data_dict
    assert MESS.SYMF in data_dict
    assert MESS.FREQ in data_dict
    assert MESS.ZENE in data_dict
    assert MESS.ELVL in data_dict
    assert MESS.CORE in data_dict
    assert MESS.MOL_PATH in path_dict
    
    # Loop over the dirs building a molecule string for each pt on the IRC
    # Build the molecule string with the core string and above data
    mol_dirs = os.listdir(path_dict[MESS.MOL_PATH])
    mol_dirs.sort(key=int)
    full_mol_dirs = [ os.path.join(path_dict[MESS.MOL_PATH], mol_dir) for mol_dir in mol_dirs]
    irc_pt_strs = []
    for i, mol_dir in enumerate(full_mol_dirs):
    
        # initialize empty string
        irc_pt_str = ''
        
        # Build the path dicts
        path_dict.update({MESS.GEOM_PATH: mol_dir+'/geom.xyz'})
        path_dict.update({MESS.FREQ_PATH: mol_dir+'/mol.harmfreq'})
        path_dict.update({MESS.MOL_ENE_PATH: mol_dir+'/mol.ene'})
        path_dict.update({MESS.MOL_ZPVE_PATH: mol_dir+'/mol.harmzpve'})
        
        # Build the molecule string using the new path dictionary and add to big str
        if i+1 == 11:
            irc_pt_str += '! IRC PT '+str(i+1)+' (SADDLEPOINT)\n'
        else:
            irc_pt_str += '! IRC PT '+str(i+1)+'\n'
        irc_pt_str += build_mess_molecule_str(head_dict, data_dict, path_dict)

        # Append to list
        irc_pt_strs.append(irc_pt_str)

    # Build the species string with the label and molecule string
    trans_state_str = mess_io.writer.write_ts_irc(head_dict[MESS.HD_LBL],
                                                  head_dict[MESS.R_LBL],
                                                  head_dict[MESS.P_LBL],
                                                  irc_pt_strs)

    return trans_state_str

"""
uses dictionaries defined in other file to build the MESS strings
"""

import mess_io
from mess_io.writer_lvl1.params import MESS
from mess_io.writer_lvl1.molecule_ import build_mess_molecule_str


def build_mess_ts_sadpt_str(trans_states):
    """ Builds the MESS input strings for each of the species
    """
    
    # Initialize the species string with the species header
    trans_states_str = mess_io.writer.stringslib.TS_HEAD_STR

    # Loop over each transition state
    for trans_state in trans_states:
    
        # Set dictionaries for clarity
        data_dict = trans_state[0]
        path_dict = trans_state[1]
        
        # SET THE VALUES FOR REQUIRED INFORMATION #
    
        # Check for existance of required keys
        assert MESS.TS_LBL in data_dict
        assert MESS.R_LBL in data_dict
        assert MESS.P_LBL in data_dict
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
        molecule_str = build_mess_molecule_str(data_dict, path_dict)
    
        # Build the species string with the label and molecule string
        trans_state_str = mess_io.writer.write_ts_sadpt(data_dict[MESS.TS_LBL],
                                                        data_dict[MESS.R_LBL],
                                                        data_dict[MESS.P_LBL],
                                                        molecule_str)

        # Add to the species string the specie string and a seperator
        trans_states_str += trans_state_str
        trans_states_str += mess_io.writer.stringslib.SPECIES_SEC_SEP_STR

    return trans_states_str
